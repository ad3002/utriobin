#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>

// Function to read the file and build a map of kmers and their frequencies
std::unordered_map<std::string, int> readKmerFile(const std::string& filename) {
    std::unordered_map<std::string, int> kmerMap;
    std::ifstream file(filename);
    std::string line;
    int lineCount = 0;

    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return kmerMap;
    }

    while (getline(file, line)) {
        std::istringstream iss(line);
        std::string kmer;
        int freq;
        if (!(iss >> kmer >> freq)) {
            continue; // error in line format, skip this line
        }
        kmerMap[kmer] = freq;
        lineCount++;
        if (lineCount % 100000000 == 0) {
            std::cout << "Read " << lineCount << " entries from " << filename << std::endl;
        }
    }
    std::cout << "Finished reading " << lineCount << " entries from " << filename << std::endl;
    return kmerMap;
}

// Main function
int main(int argc, char* argv[]) {
    if (argc != 4) {
        std::cerr << "Usage: " << argv[0] << " <file1> <file2> <output_file>" << std::endl;
        return 1;
    }

    std::string file1 = argv[1];
    std::string file2 = argv[2];
    std::string outputFile = argv[3];

    // Read the first file to build a map of kmers and their frequencies
    std::cout << "Starting to read file: " << file1 << std::endl;
    auto kmerMap1 = readKmerFile(file1);

    std::ofstream outFile(outputFile);
    if (!outFile.is_open()) {
        std::cerr << "Failed to create output file: " << outputFile << std::endl;
        return 2;
    }

    std::cout << "Processing kmers from " << file2 << " and matching against " << file1 << "..." << std::endl;
    std::ifstream file2Stream(file2);
    std::string line;
    int matchedKmers = 0;

    while (getline(file2Stream, line)) {
        std::istringstream iss(line);
        std::string kmer;
        int freq2;
        if (!(iss >> kmer >> freq2)) {
            continue; // error in line format, skip this line
        }
        auto it = kmerMap1.find(kmer);
        if (it != kmerMap1.end()) {
            int freq1 = it->second;
            outFile << kmer << "\t" << freq1 << "\t" << freq2 << std::endl;
            matchedKmers++;
            if (matchedKmers % 100000000 == 0) {
                std::cout << "Written " << matchedKmers << " matching kmers to " << outputFile << "." << std::endl;
            }
        }
    }

    outFile.close();
    std::cout << "Total matched kmers: " << matchedKmers << " (results written to " << outputFile << ")" << std::endl;
    return 0;
}

