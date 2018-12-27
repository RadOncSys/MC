#include "pch.h"
#include "DataSet.h"
#include "DicomDictionary.h"
#include "Element.h"
#include "ElementFactory.h"
#include "Tag.h"
#include <algorithm>

using namespace RO::Dicom;

#define DICOM_TRACE 0

DataSet::DataSet() 
{
	elements_.push_back(std::make_shared<std::map<UINT32, std::shared_ptr<Element>>>());
}
DataSet::DataSet(TransferSyntaxType ts) 
{ 
	elements_.push_back(std::make_shared<std::map<UINT32, std::shared_ptr<Element>>>());
	ts_ = ts; 
}

void DataSet::SetTransferSyntax(std::shared_ptr<TransferSyntaxItem> ts)
{
	ts_ = ts;
}

std::shared_ptr<Element> DataSet::GetElement(UINT32 tag, unsigned int idx)
{
	auto it = elements_[idx]->find(tag);
	if (it == elements_[idx]->end())
		return nullptr;
	else
		return it->second;
}

std::string DataSet::GetElementString(UINT32 tag, unsigned int idx)
{
	auto it = elements_[idx]->find(tag);
	if(it == elements_[idx]->end())
		return std::string();
	else
		return it->second->GetString();
}

void DataSet::AddElement(std::shared_ptr<Element> element, unsigned int idx)
{
	(*elements_[idx])[element->tag()] = element;
}

void DataSet::DeleteElement(UINT32 tag, unsigned int idx)
{
	auto it = elements_[idx]->find(tag);
	if(it != elements_[idx]->end())
		elements_[idx]->erase(it);
}

std::shared_ptr<std::vector<std::shared_ptr<Element>>> DataSet::GetElementList(unsigned int idx)
{
	auto elements = std::make_shared<std::vector<std::shared_ptr<Element>>>();
	std::for_each(elements_[idx]->begin(), elements_[idx]->end(),
		[elements](std::pair<UINT32, std::shared_ptr<Element>> x)
	{
		elements->push_back(x.second);
	});
	return elements;
}

void DataSet::CopyElementToDataSet(DataSet& ds, UINT32 tag, bool required, const char* errmsg, unsigned int idx)
{
	std::shared_ptr<Element> element = this->GetElement(tag);
	if(element)
		ds.AddElement(element, idx);
	else if(required)
		throw std::exception(errmsg);
}

void DataSet::AddElementString(UINT32 tag, const char* str, unsigned int idx)
{
	auto ditem = DicomDictionary::Item(tag);
	auto v = ElementFactory::MakeElement(tag, ditem->Vr);
	v->SetString(std::string(str));
	this->AddElement(v, idx);
}

void DataSet::AddElementInt(UINT32 tag, int i, unsigned int idx)
{
	AddElementInt(tag, std::vector<int>(1, i), idx);
}

void DataSet::AddElementInt(UINT32 tag, std::vector<int> vi, unsigned int idx)
{
	auto ditem = DicomDictionary::Item(tag);
	auto v = ElementFactory::MakeElement(tag, ditem->Vr);
	v->SetInt(vi);
	this->AddElement(v, idx);
}

void DataSet::AddElementDouble(UINT32 tag, double f, unsigned int idx)
{
	AddElementDouble(tag, std::vector<double>(1, f), idx);
}

void DataSet::AddElementDouble(UINT32 tag, std::vector<double> vf, unsigned int idx)
{
	auto ditem = DicomDictionary::Item(tag);
	auto v = ElementFactory::MakeElement(tag, ditem->Vr);
	v->SetDouble(vf);
	this->AddElement(v, idx);
}

void DataSet::AddElementDate(UINT32 tag, UINT16 year, UINT16 month, UINT16 day, unsigned int idx)
{
	auto ditem = DicomDictionary::Item(tag);
	auto v = ElementFactory::MakeElement(tag, ditem->Vr);
	v->SetDate(year, month, day);
	this->AddElement(v, idx);
}

void DataSet::AddElementDate(UINT32 tag, INT64 d, unsigned int idx)
{
	SYSTEMTIME t;
	FileTimeToSystemTime((FILETIME*)(&d), &t);
	AddElementDate(tag, t.wYear, t.wMonth, t.wDay, idx);
}

void DataSet::AddElementTime(UINT32 tag, UINT16 hour, UINT16 min, UINT16 sec, UINT32 tick, unsigned int idx)
{
	auto ditem = DicomDictionary::Item(tag);
	auto v = ElementFactory::MakeElement(tag, ditem->Vr);
	v->SetTime(hour, min, sec, tick);
	this->AddElement(v, idx);
}

void DataSet::AddElementTime(UINT32 tag, INT64 time, unsigned int idx)
{
	SYSTEMTIME t;
	FileTimeToSystemTime((FILETIME*) (&time), &t);
	AddElementTime(tag, t.wHour, t.wMinute, t.wSecond, t.wMilliseconds, idx);
}

void DataSet::AddElementDateTime(UINT32 tag, UINT16 year, UINT16 month, UINT16 day, UINT16 hour, UINT16 min, UINT16 sec, UINT32 tick, unsigned int idx)
{
	auto ditem = DicomDictionary::Item(tag);
	auto v = ElementFactory::MakeElement(tag, ditem->Vr);
	v->SetDateTime(year, month, day, hour, min, sec, tick);
	this->AddElement(v, idx);
}

void DataSet::AddElementDateTime(UINT32 tag, INT64 dt, unsigned int idx)
{
	SYSTEMTIME t;
	FileTimeToSystemTime((FILETIME*) (&dt), &t);
	AddElementDateTime(tag, t.wYear, t.wMonth, t.wDay, t.wHour, t.wMinute, t.wSecond, t.wMilliseconds, idx);
}

UINT32 DataSet::ReadTag(std::istream& is)
{
	UINT16 group = ReadUINT16(is);
	UINT16 element = ReadUINT16(is);
	return MakeTag(group, element);
}

Vr::VR DataSet::ReadVr(std::istream& is)
{
	UINT16 v = ReadUINT16(is);
	return static_cast<Vr::VR>(v);
}

UINT16 DataSet::ReadLength(std::istream& is)
{
	return ReadUINT16(is);
}

void DataSet::Read(std::istream& is, unsigned int maxlen)
{
	const char* errmsg = "error in DataSet parsing from stream";
	unsigned int count = 0;
	unsigned int idx = 0;	// индекс текущей группы для вставления элементов

	while(true)
	{
		// Check the end of data set.
		// The end could be the end of sequence or ...
		UINT32 nextValue = Element::PeekUINT32(is);
		if(nextValue == ItemDelimitationItem)
		{
			// Skip delimeters and prepare for the next dataset.
			ReadUINT32LE(is);	// delimeter
			ReadUINT32LE(is);	// must be 0

			nextValue = ReadUINT32LE(is);
			if(nextValue == SequenceDelimitationItem)
			{
				nextValue = ReadUINT32LE(is);	// must be 0
				if(nextValue != 0) 
					throw std::exception(errmsg);
				break;
			}
			else if(nextValue == Item) // begin new item of undefined length
			{
				nextValue = ReadUINT32LE(is);
				if(nextValue != 0xFFFFFFFF) 
					throw std::exception(errmsg);
			}
			else
				throw std::exception(errmsg);

			count += 16;
		}
		// В файлах GE и, кажется, согласно стандарту у итема нет флага окончания.
		// Поэтому логика выше не срабатывает и нужно повторить проверку начала нового итема.
		else if (nextValue == Item)
		{
			// HACK!! TODO: Есть принципиальная проблема.
			// Пока вся SQ читается как DataSet.
			// Но на самом деле может быть массив групп элементов.
			// А в DataSet элементы должны быть уникальными.
			// А группа не имеет своего тега.
			// Значит модель DICOM файла должна быть расширена на поддержку таких групп.

			// Временно просто пропускаем разделители чтобы чтение не вываливалось.
			ReadUINT32(is);		// UINT32 itemDelimeter
			ReadUINT32(is);		// UINT32 itemlen

			// Это место, где появляется новая группа элеменотов (не DataSet).
			// Добавляем группу и обновляем индекс, куда нужно помещать поступающие элементы.
			// Это решение проблемы выше.
			elements_.push_back(std::make_shared<std::map<UINT32, std::shared_ptr<Element>>>());
			idx = (unsigned	int)elements_.size() - 1;
#if DICOM_TRACE > 0
			wchar_t txt[64];
			swprintf_s(txt, 64, L"New element group eded, idx = %i\n", idx);
			OutputDebugString(txt);
#endif
			count += 8;
		}

		UINT32 tag = ReadTag(is);
		if(is.fail()) throw std::exception(errmsg);
		count += 4;

		UINT32 len;

		Vr::VR vr;
		if(ts_->ExplicitVr)
		{
			vr = ReadVr(is);
			if(is.fail()) throw std::exception(errmsg);
			len = ReadUINT16(is);
			if(is.fail()) throw std::exception(errmsg);
			count += 4;
		}
		else
		{
			auto ditem = DicomDictionary::Item(tag);
			if(!ditem)
				vr = Vr::SH;	// HACK !!!
			else
				vr = ditem->Vr;
			len = ReadUINT32(is);
			if(is.fail()) throw std::exception(errmsg);
			count += 4;
		}

		auto v = ElementFactory::MakeElement(tag, vr);

#if DICOM_TRACE > 0
		if (vr == Vr::SQ)
			OutputDebugString((std::wstring(L"SQ \t") + ((Element*)v.get())->ToWString() + L"\n").c_str());
#endif

		v->Read(is, len, ts_);
		count += len;

		AddElement(v, idx);

#if DICOM_TRACE > 0
		if (vr == Vr::SQ)
			OutputDebugString(L"Eof SQ\n");
		else
			OutputDebugString((v->ToWString() + L"\n").c_str());
#endif
		
		if (maxlen != 0xFFFFFFFF && count >= maxlen)
		{
			// Обработка окончания объекта при ограничении по длине.
			// Реально имеет отношение только к SQ.
			// Означает окончание ценпочки элементов, когда не используется ограничитель
			break;
		}

		// Pixel data should be the last element without any delimeters
		if(tag == PixelData)
			break;
	}
}

void DataSet::WriteTag(std::ostream& os, UINT32 tag)
{
	WriteUINT16(os, TagGroup(tag));
	WriteUINT16(os, TagElement(tag));
}

void DataSet::WriteVr(std::ostream& os, Vr::VR vr)
{
	WriteUINT16(os, static_cast<UINT16>(vr));
}

void DataSet::WriteLength(std::ostream& os, UINT16 len)
{
	WriteUINT16(os, len);
}

void DataSet::Write(std::ostream& os, TransferSyntaxType ts)
{
	ts_ = ts;

	for each (auto elements in elements_)
	{
		for_each(elements->begin(), elements->end(),
			[this, &os](std::pair<UINT32, std::shared_ptr<Element>> x)
		{
			WriteTag(os, x.first);
			if (ts_->ExplicitVr)
				WriteVr(os, x.second->vr());
			x.second->Write(os, ts_);

			if (os.fail())
				throw std::exception("error in DataSet writing to stream");
		});
	}
}
