#include "ParseHandler.h"


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


DataHolder PipeReader()
{
	
	LOG(INFO) << "Detecting Data Stream. Switching to Piped Input Mode";
	if (Settings.DataFile != "_no_file_")
	{
		LOG(WARN) << "Multiple input methods provided (input stream & '" << Settings.DataFile << "').\n\tData Stream has priority";
	}
	std::ios_base::sync_with_stdio(false);
	std::cin.rdbuf()->pubsetbuf(new char[1024 * 1024], 1024 * 1024); // 1 MB buffer
	return ParseRawInput(std::cin);

}

DataHolder ShellExecute(std::string cmd)
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
	auto out = ParseRawInput(stream);
	auto exit = pclose(pipe);
	if (WEXITSTATUS(exit) != 0)
	{
		throw std::runtime_error("Command (" + cmd +") returned a non-zero exit code");
	}
	return out;
}


DataHolder ArchiveReader()
{
	JAR::Archive archive(Settings.DataFile,std::ios::in);

	auto archiveFiles = archive.ListFiles();
	if (archiveFiles.size() == 0)
	{
		throw runtime_error("There are no files found within the requested datafile '" + Settings.DataFile + "'.\n\tAre you sure it is a valid GenCHORD Archive file?");
	}
	std::ostringstream log;
	log << "Archive found " << archiveFiles.size() << " files within";
	bool manifestPresent = false;
	for (auto file : archiveFiles)
	{
		if (file == MANIFEST_FILE_NAME)
		{
			manifestPresent = true;
		}
		log << "\n\t" << file;
	}
	LOG(DEBUG) << log.str();

	std::vector<std::string> relevantFiles;
	if (manifestPresent)
	{
		LOG(DEBUG) << "Found manifest file";
		
		std::vector<std::string> chroms = JSL::split(archive.Text(MANIFEST_FILE_NAME),'\n');
		for (auto q : chroms)
		{
			std::string targetFile = q + "_" + std::to_string(Settings.AccumulationFactor) + ".dat";
			if (JSL::FindXInY(targetFile,archiveFiles) == -1)
			{
				throw runtime_error(Settings.DataFile + " does not contain the file " + targetFile + "\n\tThere are two potential causes for this.\n\t\t1. A corrupted or incomplete archive \n\t\t2. A non-archived accumulation factor (" + std::to_string(Settings.AccumulationFactor) + ").\n\tYou will need to rerun with the original .bam/.sam file");
			}

			relevantFiles.push_back(targetFile);
		}
	}
	else
	{
		LOG(WARN)  << "WARNING! No manifest file (" << MANIFEST_FILE_NAME << ") found within the archive " << Settings.DataFile << "\n\tChromosomes will be analysed in a random order";

		std::regex re("[_.]");
		for (auto file : archiveFiles)
		{
			
			std::sregex_token_iterator first{file.begin(), file.end(), re, -1}, last;//the '-1' is what makes the regex split (-1 := what was not matched)
			std::vector<std::string> tokens{first, last};
			if (tokens.size() == 3 && std::stoi(tokens[1]) == Settings.AccumulationFactor)
			{
				relevantFiles.push_back(file);
			}
		}
	}

	if (relevantFiles.size() == 0)
	{
		throw runtime_error("No files were found in the archive which matched your specifications");
	}

	DataHolder output;
	std::vector<std::tuple<dnaindex,lint, lint>> chromVector;
	LOG(DEBUG) << "Beginning file read";
	for (auto file: relevantFiles)
	{
		LOG(DEBUG) << "\tOpening file " << file;
		auto name = JSL::split(file,'_')[0];
		archive.ReadTabular(file,chromVector,' ');
		output.Append(name,chromVector);
	}
	return output;
}

DataHolder FileReader()
{
	LOG(DEBUG) << "File input detected, determining how to open file";

	if (!StringIsSanitised(Settings.DataFile))
	{
		throw std::runtime_error("The datafile \'" + Settings.DataFile + "\' contains unsafe characters (|;& etc), and therefore cannot be opened.");
	}
	auto file = Settings.DataFile;
	
	std::string fileExtension = JSL::split(file,'.').back();
	if (fileExtension == "bam" || fileExtension == "cram" || fileExtension == "sam")
	{
		LOG(INFO) << "." << fileExtension << " extension detected, calling samtools depth";
		return ShellExecute("samtools depth " + file);
	}
	else if (fileExtension == "gca")
	{
		LOG(INFO) << ".gca (GenCHORD Archive) file detected. Extracting";
		if (Settings.ProcessMode)
		{
			LOG(WARN) << ".gca files are already archives. PROCESS MODE has no effect on existing archives";
		}
		return ArchiveReader();
	}
	else
	{
		LOG(WARN) << "Unrecognised file extension. Attempting to read as text file";
		return ShellExecute("cat " + file);
	}
}


DataHolder ParseData()
{

	LOG(INFO) << "Beginning parsing of input data";

	if (JSL::PipedInputFound())
	{
		return PipeReader();
	}
	else if (Settings.DataFile != "_no_file_")
	{
		return FileReader();
	}
	else
	{
		throw runtime_error("No data was provided for analysis.");
	}

}