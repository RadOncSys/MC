// THIS CODE BELONGS TO Radiation Oncology Intellectual Systems and Services LLC
// Copyright (c) 2013, ROISS LLC. All rights reserved
//
// Author: Gennady Gorlachev (ggorlachev@roiss.ru)
//---------------------------------------------------------------------------
#pragma once
#include "RO.Dicom.h"
#include "Dataset.h"
#include "TransferSyntax.h"
#include <memory>
#include <vector>
#include <string>

namespace RO
{
	namespace Dicom
	{
		class Dataset;

		RODICOM_API UINT16 ReadUINT16(std::istream& is);
		RODICOM_API UINT32 ReadUINT32(std::istream& is);
		RODICOM_API UINT32 ReadUINT32LE(std::istream& is);

		RODICOM_API void WriteUINT16(std::ostream& os, UINT16 v);
		RODICOM_API void WriteUINT32(std::ostream& os, UINT32 v);
		RODICOM_API void WriteUINT32LE(std::ostream& os, UINT32 v);

		// Element Base Class
		class RODICOM_API Element
		{
		public:
			Element(UINT32 tag, enum Vr::VR vr) : tag_(tag), vr_(vr) {}
			virtual ~Element() {}

			UINT32 tag() const { return tag_; }
			enum Vr::VR vr() const { return vr_; }

			// Each Element type class should have its own read and write functions
			// to get and set data in stream
			virtual void Read(std::istream& is, UINT32 len, TransferSyntaxType ts) = 0;
			virtual void Write(std::ostream& os, TransferSyntaxType ts) const = 0;

			// Base class will normaly throw exception on the following functions.
			// Those Element classes, which can manage given data types
			// will owerride corresponding virtual functions to provide neccessary functionality.
			virtual std::string GetString() const;
			virtual void SetString(const std::string& s);

			virtual std::vector<double> GetDouble() const;
			virtual void SetDouble(std::vector<double>);

			virtual std::vector<int> GetInt() const;
			virtual void SetInt(std::vector<int>);

			virtual std::shared_ptr<DataSet> GetDataSet() const;
			virtual void SetDataSet(std::shared_ptr<DataSet>);

			// Date and time functions
			virtual void GetDate(UINT16& year, UINT16& month, UINT16& day) const;
			virtual void GetTime(UINT16& hour, UINT16& min, UINT16& sec, UINT32& tick) const;
			virtual void GetDateTime(UINT16& year, UINT16& month, UINT16& day, UINT16& hour, UINT16& min, UINT16& sec, UINT32& tick) const;
			virtual void SetDate(UINT16 year, UINT16 month, UINT16 day);
			virtual void SetTime(UINT16 hour, UINT16 min, UINT16 sec, UINT32 tick);
			virtual void SetDateTime(UINT16 year, UINT16 month, UINT16 day, UINT16 hour, UINT16 min, UINT16 sec, UINT32 tick);

			// DEBUG
			virtual std::wstring ToWString() const;

			static UINT32 PeekUINT32(std::istream& is);
			static UINT32 ReadDelimeter(std::istream& is);

		protected:
			// Common for all element types utilities
			static void ReadTag(std::istream& is);
			static void ReadVr(std::istream& is);
			static long ReadLength(std::istream& is);

			static void WriteTag(std::ostream& os);
			static void WriteVr(std::ostream& os);
			static void WriteLength(std::ostream& os, long lengts);

		private:
			UINT32 tag_;
			enum Vr::VR vr_;
		};

		// String
		class RODICOM_API ElementString : public Element
		{
		public:
			ElementString(UINT32 tag, enum Vr::VR vr) : Element(tag, vr) {}

			void Read(std::istream& is, UINT32 len, TransferSyntaxType ts) override;
			void Write(std::ostream& os, TransferSyntaxType ts) const override;
			std::string GetString() const override;
			void SetString(const std::string& s) override;
			virtual std::vector<int> GetInt() const override;
			virtual void SetInt(std::vector<int>) override;
			std::vector<double> GetDouble() const override;
			void SetDouble(std::vector<double>) override;
			virtual std::wstring ToWString() const;
		protected:
			std::string data_;
		};

		// OtherString
		class RODICOM_API ElementOtherString : public Element
		{
		public:
			ElementOtherString(UINT32 tag, enum Vr::VR vr) : Element(tag, vr), data_(nullptr), len_(0) {}
			virtual ~ElementOtherString();
			void Read(std::istream& is, UINT32 len, TransferSyntaxType ts) override;
			void Write(std::ostream& os, TransferSyntaxType ts) const override;

			std::string GetString() const override;
			void SetString(const std::string& s) override;
			std::vector<int> GetInt() const override;
			void SetInt(std::vector<int>) override;
			std::vector<double> GetDouble() const override;
			void SetDouble(std::vector<double>) override;

			const char* GetData() { return data_; }
			int GetLength() const { return len_; }

			void SetData(const char* data, int len);
			virtual std::wstring ToWString() const;

		protected:
			char* data_;
			int len_;
		};

		// AE
		class RODICOM_API ElementAE : public ElementString
		{
		public:
			ElementAE(UINT32 tag, enum Vr::VR vr) : ElementString(tag, vr) {}
			void SetString(const std::string& s) override;
		};

		// AS
		class RODICOM_API ElementAS : public ElementString
		{
		public:
			ElementAS(UINT32 tag, enum Vr::VR vr) : ElementString(tag, vr) {}
			void SetString(const std::string& s) override;
		};

		// AT
		class RODICOM_API ElementAT : public Element
		{
		public:
			ElementAT(UINT32 tag, enum Vr::VR vr) : Element(tag, vr) {}

			void Read(std::istream& is, UINT32 len, TransferSyntaxType ts) override;
			void Write(std::ostream& os, TransferSyntaxType ts) const override;

			std::vector<int> GetInt() const override;
			void SetInt(std::vector<int>) override;
			virtual std::wstring ToWString() const;

		private:
			std::vector<int> data_;
		};

		// CS
		class RODICOM_API ElementCS : public ElementString
		{
		public:
			ElementCS(UINT32 tag, enum Vr::VR vr) : ElementString(tag, vr) {}
			void SetString(const std::string& s) override;
		};

		// DA
		class RODICOM_API ElementDA : public ElementString
		{
		public:
			ElementDA(UINT32 tag, enum Vr::VR vr) : ElementString(tag, vr) {}
			void GetDate(UINT16& year, UINT16& month, UINT16& day) const override;
			void SetDate(UINT16 year, UINT16 month, UINT16 day) override;
		};

		// TM
		class RODICOM_API ElementTM : public ElementString
		{
		public:
			ElementTM(UINT32 tag, enum Vr::VR vr) : ElementString(tag, vr) {}
			void GetTime(UINT16& hour, UINT16& min, UINT16& sec, UINT32& tick) const override;
			void SetTime(UINT16 hour, UINT16 min, UINT16 sec, UINT32 tick) override;
		};

		// DT
		class RODICOM_API ElementDT : public ElementString
		{
		public:
			ElementDT(UINT32 tag, enum Vr::VR vr) : ElementString(tag, vr) {}
			void GetDateTime(UINT16& year, UINT16& month, UINT16& day, UINT16& hour, UINT16& min, UINT16& sec, UINT32& tick) const override;
			void SetDateTime(UINT16 year, UINT16 month, UINT16 day, UINT16 hour, UINT16 min, UINT16 sec, UINT32 tick) override;
		};

		// FL
		class RODICOM_API ElementFL : public Element
		{
		public:
			ElementFL(UINT32 tag, enum Vr::VR vr) : Element(tag, vr) {}
			void Read(std::istream& is, UINT32 len, TransferSyntaxType ts) override;
			void Write(std::ostream& os, TransferSyntaxType ts) const override;
			std::vector<double> GetDouble() const override;
			void SetDouble(std::vector<double>) override;
			virtual std::wstring ToWString() const;
		private:
			std::vector<double> data_;
		};

		// FD
		class RODICOM_API ElementFD : public Element
		{
		public:
			ElementFD(UINT32 tag, enum Vr::VR vr) : Element(tag, vr) {}
			void Read(std::istream& is, UINT32 len, TransferSyntaxType ts) override;
			void Write(std::ostream& os, TransferSyntaxType ts) const override;
			std::vector<double> GetDouble() const override;
			void SetDouble(std::vector<double>) override;
			virtual std::wstring ToWString() const;
		private:
			std::vector<double> data_;
		};

		// SL
		class RODICOM_API ElementSL : public Element
		{
		public:
			ElementSL(UINT32 tag, enum Vr::VR vr) : Element(tag, vr) {}
			void Read(std::istream& is, UINT32 len, TransferSyntaxType ts) override;
			void Write(std::ostream& os, TransferSyntaxType ts) const override;
			std::vector<int> GetInt() const override;
			void SetInt(std::vector<int>) override;
			virtual std::wstring ToWString() const;
		private:
			std::vector<int> data_;
		};

		// SS
		class RODICOM_API ElementSS : public Element
		{
		public:
			ElementSS(UINT32 tag, enum Vr::VR vr) : Element(tag, vr) {}
			void Read(std::istream& is, UINT32 len, TransferSyntaxType ts) override;
			void Write(std::ostream& os, TransferSyntaxType ts) const override;
			std::vector<int> GetInt() const override;
			void SetInt(std::vector<int>) override;
			virtual std::wstring ToWString() const;
		private:
			std::vector<int> data_;
		};

		// UI
		class RODICOM_API ElementUI : public ElementString
		{
		public:
			ElementUI(UINT32 tag, enum Vr::VR vr) : ElementString(tag, vr) {}
			void SetString(const std::string& s) override;
		};

		// UL
		class RODICOM_API ElementUL : public Element
		{
		public:
			ElementUL(UINT32 tag, enum Vr::VR vr) : Element(tag, vr) {}
			void Read(std::istream& is, UINT32 len, TransferSyntaxType ts) override;
			void Write(std::ostream& os, TransferSyntaxType ts) const override;
			std::vector<int> GetInt() const override;
			void SetInt(std::vector<int>) override;
			virtual std::wstring ToWString() const;
		private:
			std::vector<int> data_;
		};

		// US
		class RODICOM_API ElementUS : public Element
		{
		public:
			ElementUS(UINT32 tag, enum Vr::VR vr) : Element(tag, vr) {}

			void Read(std::istream& is, UINT32 len, TransferSyntaxType ts) override;
			void Write(std::ostream& os, TransferSyntaxType ts) const override;

			std::vector<int> GetInt() const override;
			void SetInt(std::vector<int>) override;
			virtual std::wstring ToWString() const;

		private:
			std::vector<int> data_;
		};

		// SQ
		class RODICOM_API ElementSQ : public Element
		{
		public:
			ElementSQ(UINT32 tag, enum Vr::VR vr) : Element(tag, vr) {}
			void Read(std::istream& is, UINT32 len, TransferSyntaxType ts) override;
			void Write(std::ostream& os, TransferSyntaxType ts) const override;
			std::shared_ptr<DataSet> GetDataSet() const override;
			void SetDataSet(std::shared_ptr<DataSet>) override;
		private:
			std::shared_ptr<DataSet> data_;
		};

	}
}
