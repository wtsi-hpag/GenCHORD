#include "Streamer.h"
#include "Data.h"
#include "Aggregators.h"



// Helper functor to adapt `FILE*` to `std::istream`
class PopenBuffer : public std::streambuf
{
	public:
		PopenBuffer(FILE* file) : file_(file) {}
	protected:
		int underflow() override
		{
			if (gptr() == egptr())
			{
				size_t n = fread(buffer_, 1, sizeof(buffer_), file_);
				if (n == 0) return traits_type::eof();
				setg(buffer_, buffer_, buffer_ + n);
			}
			return traits_type::to_int_type(*gptr());
		}
	private:
		FILE* file_;
		char buffer_[1024];
};


void PipeReader()
{
	LOG(INFO) << "Detecting Data Stream. Switching to Piped Input Mode";
	if (Settings.DataFile != "_no_file_")
	{
		LOG(WARN) << "Multiple input methods provided (input stream & '" << Settings.DataFile << "').\n\tData Stream has priority";
	}

	AggregateStream(std::cin);

}

void ShellExecute(std::string cmd)
{
	// std::string cmd = "cat Data/Aaron.dat";
	LOG(DEBUG) << "Calling popen with command '" << cmd << "'";
	FILE* pipe = popen(cmd.c_str(),"r");
	if (!pipe)
	{
		throw std::runtime_error("Failed to open pipe for command: " + cmd);
	}
	PopenBuffer tmp(pipe);
	std::istream stream(&tmp);
	AggregateStream(stream);
	auto exit = pclose(pipe);
	if (WEXITSTATUS(exit) != 0)
	{
		throw std::runtime_error("Command (" + cmd +") returned a non-zero exit code");
	}
}

void FileReader()
{
	LOG(INFO) << "File input detected, determining how to open file";

	if (!StringIsSanitised(Settings.DataFile))
	{
		throw std::runtime_error("The datafile \'" + Settings.DataFile + "\' contains unsafe characters (|;& etc), and therefore cannot be opened.");
	}
	auto file = Settings.DataFile;
	
	std::string fileExtension = JSL::split(file,'.').back();
	if (fileExtension == "bam")
	{
		LOG(INFO) << ".bam extension detected, calling samtools depth";
		ShellExecute("samtools depth " + file);
	}
	else if (fileExtension == "gca")
	{
		LOG(INFO) << ".gca (GenCHORD Archive) file detected. Extracting";
	}
	else
	{
		LOG(WARN) << "Unrecognised file extension. Attempting to read as text file";
		ShellExecute("cat " + file);
	}


	
	
}
void ParseData()
{
	LOG(INFO) << "Beginning parsing of input data";

	if (JSL::PipedInputFound())
	{
		PipeReader();
	}
	else if (Settings.DataFile != "_no_file_")
	{
		FileReader();
	}
	else
	{
		LOG(WARN) << "No data was provided for analysis.";
	}

}