#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>

// Function to read the original kmer file with two columns
std::unordered_map<std::string, int> readKmerFile(const std::string& filename) {
    std::unordered_map<std::string, int> kmerMap;
    std::ifstream file(filename);
    std::string line;
    size_t lineCount = 0;

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

// Function to read computed kmer files with three columns
std::unordered_map<std::string, std::pair<int, int>> readComputedKmerFile(const std::string& filename, const std::unordered_map<std::string, int>& filterMap) {
    std::unordered_map<std::string, std::pair<int, int>> kmerMap;
    std::ifstream file(filename);
    std::string line;
    size_t lineCount = 0;

    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return kmerMap;
    }

    while (getline(file, line)) {
        std::istringstream iss(line);
        std::string kmer;
        int freq1, freq2;
        if (!(iss >> kmer >> freq1 >> freq2)) {
            continue; // error in line format, skip this line
        }
        if (filterMap.find(kmer) != filterMap.end()) {
            kmerMap[kmer] = std::make_pair(freq1, freq2);
            lineCount++;
            if (lineCount % 100000000 == 0) {
                std::cout << "Processed " << lineCount << " relevant entries from " << filename << std::endl;
            }
        }
    }
    std::cout << "Finished processing relevant entries from " << filename << std::endl;
    return kmerMap;
}

// Main function
int main(int argc, char* argv[]) {
    if (argc != 5) {
        std::cerr << "Usage: " << argv[0] << " <original_file> <computed_file1> <computed_file2> <output_file>" << std::endl;
        return 1;
    }

    std::string originalFile = argv[1];
    std::string computedFile1 = argv[2];
    std::string computedFile2 = argv[3];
    std::string outputFile = argv[4];

    auto originalKmers = readKmerFile(originalFile);
    auto computedKmers1 = readComputedKmerFile(computedFile1, originalKmers);
    auto computedKmers2 = readComputedKmerFile(computedFile2, originalKmers);

    std::ofstream outFile(outputFile);
    if (!outFile.is_open()) {
        std::cerr << "Failed to create output file: " << outputFile << std::endl;
        return 2;
    }

    // Write headers
    outFile << "Kmer\tFreq1 (File1)\tFreq2 (File1)\tFreq1 (File2)\tFreq2 (File2)" << std::endl;

    size_t processedKmers = 0;
    for (const auto& pair : originalKmers) {
        const std::string& kmer = pair.first;
        auto it1 = computedKmers1.find(kmer);
        auto it2 = computedKmers2.find(kmer);

        if (it1 != computedKmers1.end() && it2 != computedKmers2.end()) {
            outFile << kmer << "\t"
                    << it1->second.first << "\t" << it1->second.second << "\t"
                    << it2->second.first << "\t" << it2->second.second << std::endl;
            processedKmers++;
            
        } else if (it1 != computedKmers1.end()) {
            outFile << kmer << "\t"
                    << it1->second.first << "\t" << it1->second.second << "\t"
                    << 0 << "\t" << 0 << std::endl;
            processedKmers++;
        } else if (it2 != computedKmers1.end()) {
            outFile << kmer << "\t"
                    << 0 << "\t" << 0 << "\t"
                    << it2->second.first << "\t" << it2->second.second << std::endl;
            processedKmers++;
        } else {
            outFile << kmer << "\t"
                    << 0 << "\t" << 0 << "\t"
                    << 0 << "\t" << 0 << std::endl;
            processedKmers++;

        }

        if (processedKmers % 100000000 == 0) {
            std::cout << "Written " << processedKmers << " matching kmers to " << outputFile << "." << std::endl;
        }
    }

    outFile.close();
    std::cout << "Total processed kmers: " << processedKmers << " (results written to " << outputFile << ")" << std::endl;
    return 0;
}

