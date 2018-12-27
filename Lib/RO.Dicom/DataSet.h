// THIS CODE BELONGS TO Radiation Oncology Intellectual Systems and Services LLC
// Copyright (c) 2013, ROISS LLC. All rights reserved
//
// Author: Gennady Gorlachev (ggorlachev@roiss.ru)
//---------------------------------------------------------------------------
#pragma once
#include "RO.Dicom.h"
#include "Vr.h"
#include "TransferSyntax.h"
#include <map>
#include <memory>
#include <vector>

namespace RO
{
	namespace Dicom
	{
		class Element;

		// DataSet format is described in the standard:
		//Section 7 The Data Set (PS 3.5-2011 Page 35)

		class RODICOM_API DataSet
		{
		public:
			DataSet();
			DataSet(TransferSyntaxType ts);

			void SetTransferSyntax(TransferSyntaxType ts);

			std::shared_ptr<Element> GetElement(UINT32 tag, unsigned int idx = 0);

			/// Get element string value even if element is not available.
			/// Empty string will be returned in last case.
			/// If element exist and do not support string conversion, exception will be thrown.
			std::string GetElementString(UINT32 tag, unsigned int idx = 0);

			// If element with the given tag already exista, it will be replaced.
			void AddElement(std::shared_ptr<Element>, unsigned int idx = 0);

			void DeleteElement(UINT32 tag, unsigned int idx = 0);

			unsigned int NElementGroups() const { return (unsigned int)elements_.size(); }
			std::shared_ptr<std::map<UINT32, std::shared_ptr<Element>>> Elements(unsigned int idx = 0) { return elements_[idx]; }
			std::shared_ptr<std::vector<std::shared_ptr<Element>>> GetElementList(unsigned int idx = 0);

			void CopyElementToDataSet(DataSet& ds, UINT32 tag, bool required, const char* errmsg, unsigned int idx = 0);

			// Utility methods
			void AddElementString(UINT32 tag, const char* str, unsigned int idx = 0);
			void AddElementInt(UINT32 tag, int i, unsigned int idx = 0);
			void AddElementInt(UINT32 tag, std::vector<int> vi, unsigned int idx = 0);
			void AddElementDouble(UINT32 tag, double f, unsigned int idx = 0);
			void AddElementDouble(UINT32 tag, std::vector<double> vf, unsigned int idx = 0);
			void AddElementDate(UINT32 tag, UINT16 year, UINT16 month, UINT16 day, unsigned int idx = 0);
			void AddElementDate(UINT32 tag, INT64 date, unsigned int idx = 0);
			void AddElementTime(UINT32 tag, UINT16 hour, UINT16 min, UINT16 sec, UINT32 tick, unsigned int idx = 0);
			void AddElementTime(UINT32 tag, INT64 time, unsigned int idx = 0);
			void AddElementDateTime(UINT32 tag, UINT16 year, UINT16 month, UINT16 day, UINT16 hour, UINT16 min, UINT16 sec, UINT32 tick, unsigned int idx = 0);
			void AddElementDateTime(UINT32 tag, INT64 time, unsigned int idx = 0);

			// Complete dataset (sequence is a dataset).
			// Sequences are readed reqursively
			// Если maxlen = 0, то длина неизвестна и нужно отслеживать метки окончания элементов.
			// Проблема в том, что если длина указана, то меток ограничений нет.
			void Read(std::istream& is, unsigned int maxlen = 0xFFFFFFFF);
			void Write(std::ostream& os, TransferSyntaxType ts);

			static UINT32 ReadTag(std::istream& is);
			static Vr::VR ReadVr(std::istream& is);
			static UINT16 ReadLength(std::istream& is);

			static void WriteTag(std::ostream& os, UINT32 tag);
			static void WriteVr(std::ostream& os, Vr::VR vr);
			static void WriteLength(std::ostream& os, UINT16 len);

		private:
			std::shared_ptr<TransferSyntaxItem> ts_;
			std::vector<std::shared_ptr<std::map<UINT32, std::shared_ptr<Element>>>> elements_;
		};
	}
}
