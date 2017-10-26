
#include "AADnumber.h"

const struct Number::LeafType Number::leaf;
const struct Number::UnaryType Number::unary;
const struct Number::BinaryType Number::binary;

Tape globalTape;
thread_local Tape* Number::tape = &globalTape;