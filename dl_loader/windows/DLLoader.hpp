#pragma once
#ifndef IDL_LOADER_WIN_HPP
#define IDL_LOADER_WIN_HPP

#ifdef WIN32

#include <Windows.h>

#include "../IDLLoader.hpp"

namespace dlloader {

template <typename T>
class DLLoader : public IDLLoader<T> {
  HMODULE handle;
  std::string pathToLib;
  std::string allocClassSymbol;
  std::string deallocClassSymbol;

 public:
  DLLoader(const std::string &pathToLib, const std::string &allocClassSymbol,
           const std::string &deallocClassSymbol)
      : handle(nullptr),
        pathToLib(pathToLib),
        allocClassSymbol(allocClassSymbol),
        deallocClassSymbol(deallocClassSymbol) {}

  ~DLLoader() = default;

  void DLOpenLib() override {
    handle = LoadLibrary(pathToLib.c_str());

    if (!handle) {
      std::cerr << "Error opening library: " << pathToLib << std::endl;
    }
  }

  std::shared_ptr<T> DLGetInstance() override {
    using allocClass = T *(*)();
    using deallocClass = void (*)(T *);

    auto allocFunc = reinterpret_cast<allocClass>(
        GetProcAddress(handle, allocClassSymbol.c_str()));

    auto deleteFunc = reinterpret_cast<deallocClass>(
        GetProcAddress(handle, deallocClassSymbol.c_str()));

    if (!allocFunc || !deleteFunc) {
      DLCloseLib();
      std::cerr << "Error loading symbols: " << pathToLib << std::endl;
    }

    return std::shared_ptr<T>(allocFunc(),
                              [deleteFunc](T *ptr) { deleteFunc(ptr); });
  }

  void DLCloseLib() override {
    if (FreeLibrary(handle) == 0) {
      std::cerr << "Error closing library: " << pathToLib << std::endl;
    }
  }
};
}  // namespace dlloader

#endif  // _WIN32

#endif  // IDL_LOADER_WIN_HPP