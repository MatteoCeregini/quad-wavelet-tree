#include <iostream>
#include <string>
#include <filesystem>
#include "essentials.hpp"


int main(int argc, char** argv) {
    if (argc != 3) {
        std::cout << "USAGE: ./create_dataset PATH_TO_FILE PATH_TO_DATASET" << std::endl;
        return -1;
    }
    const char* path_to_file = argv[1];
    const char* path_to_dataset = argv[2];

    // Check if the file exists
    if (!std::filesystem::exists(path_to_file)) {
        std::cout << "ERROR: \'" << path_to_file << "\' does not exists!" << std::endl;
        return -1;
    }   
    std::cout << "[1] Creating the dataset object..." << std::endl;
    essentials::dataset<uint16_t> ds; // the file is encoded as as sequence of 16-bit integers.
    
    // Try to create the dataset from the file
    if (!essentials::create_dataset(ds, path_to_file)) {
        std::cout << "ERROR: can't open: \'" << path_to_file << "\'" << std::endl;
        return -1;
    }
    std::cout << "[2] Saving the dataset object to disk..." << std::endl;
    essentials::save(ds, path_to_dataset);
    return 0;
}