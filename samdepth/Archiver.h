// The Write command and WriteMetaData structure are (verbatim, aside from some name changes) derived from 
// the `tar_to_stream` project by Armchair-Software, available at:
// https://github.com/Armchair-Software/tar_to_stream/
//
// The original code is provided under the MIT License (https://opensource.org/licenses/MIT).
// The Read functions and object-oriented implementation were written by JFG-2025.

#pragma once
#include <span>
#include <string_view>
#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>
#include <functional>
#include <cstring>
#include <stdexcept>

namespace JAR
{
	struct ReadMetaData
	{
		std::string filename;
		size_t file_size;
		std::streampos data_offset; // Byte offset to the file's data in the archive
	};

	struct WriteMetaData {
		/// Properties of the file to enter into the stream
		std::string const &filename;                                                  /// name of the file to write
		std::span<std::byte const> data;                                              /// the location of the file's contents in memory
		uint64_t mtime{0u};                                                           /// file modification time, in seconds since epoch
		std::string filemode{"644"};                                                  /// file mode
		unsigned int uid{0u};                                                         /// file owner user ID
		unsigned int gid{0u};                                                         /// file owner group ID
		std::string const &uname{"root"};                                             /// file owner username
		std::string const &gname{"root"};                                             /// file owner group name
	};


	class Archive
	{
		private:
			std::fstream Stream;
			std::unordered_map<std::string, ReadMetaData> FileIndex;
			constexpr static size_t BLOCK_SIZE = 512;
			bool IndexBuilt;
			bool HasWritten;
			void BuildIndex();
			bool ReadBlock(char*buffer);
			void WriteCleanup(unsigned int tail_length = 512u * 2u);
			// void ActivateStream(std::string archivePath, std::ios_base::openmode mode);
		public:
			Archive(std::string archivePath);
			Archive(std::string archivePath, std::ios_base::openmode mode);

			~Archive();
			std::vector<std::string> ListFiles();
			void Write(WriteMetaData &&input);
			void Write(const std::string & fileName, const std::string & data);

			Archive(const Archive&) = delete;
			Archive& operator=(const Archive&) = delete;
	};
}