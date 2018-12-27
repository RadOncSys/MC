#include "pch.h"
#include "Uid.h"
#include <cctype>
#include <strstream>

using namespace RO::Dicom;

void StripTrailingNull(std::string& str)
{
	std::string::size_type last=str.find_last_not_of('\0');
	if(last==std::string::npos)
		last=0;
	else
		last++;
	str.resize(last);
}

void StripTrailingWhitespace(std::string& str)
{
	std::string::size_type last=str.find_last_not_of(' ');
	if(last==std::string::npos)
		last=0;
	else
		last++;
	str.resize(last);

	StripTrailingNull(str);
}

void ThrowIfInvalid(const char c)
{
	if(c!='.' && !std::isdigit(c) && c!='*' )
	{
		std::strstream s;
		s << c << " is an invalid UID character";
		throw std::exception(s.str());
	}
}

/*
Uid::Uid(const std::string& s) : data_(s)
{
	if(data_.size()>64)
	{
		std::strstream error;
		error << "UID : " << s << " is too long";
		throw std::exception(error.str());
	}
	StripTrailingWhitespace(data_);
	//StripTrailingNull(data_);
	std::for_each(data_.begin(),data_.end(),ThrowIfInvalid);

}

std::string Uid::str() const
{
	return data_;
}

Uid makeUID(const std::string &Prefix)
{
	std::strstream stream;
	stream	<< Prefix
		<< "." << time(0)
		<< "." << clock()
		<< "." << rand();
	std::string s=stream.str();
	if(s.size()>64)
		throw std::exception("Generated UID is larger than 64 characters.");
	return Uid(s);
}
*/
