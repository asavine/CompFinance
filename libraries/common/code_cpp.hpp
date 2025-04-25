// Modified by Bernardo Cohen
// Original comes from
// https://caiorss.github.io/C-Cpp-Notes/CwrapperToQtLibrary.html

#pragma once
#ifndef MY_LIBRRARY_CODE_CPP_HPP
#define MY_LIBRRARY_CODE_CPP_HPP

#ifdef __cplusplus
#define EXPORT_C extern "C"
#else
#define EXPORT_C
#endif

#endif  // MY_LIBRRARY_CODE_CPP_HPP