#ifndef _PREPROCESSOR_HPP
#define _PREPROCESSOR_HPP

#define TO_STRING_INTERNAL(x) #x

/// Convert to a string
#define TO_STRING(x) \
  TO_STRING_INTERNAL(x)

/// Concatenate two tokens
///
/// Internal implementation
#define _CONCAT2_IMPL(X,Y)			\
  X ## Y

/// Concatenate two tokens
///
///  Wrapper to beat CPP
#define _CONCAT2_WRAP(X,Y)			\
  _CONCAT2_IMPL(X,Y)

/// Concatenate two tokens
///
/// User accessible implementation
#define CONCAT2(X,Y)				\
  _CONCAT2_WRAP(X,Y)

/// Concatenate three tokens
#define CONCAT3(X,Y,Z)				\
  CONCAT2(CONCAT2(X,Y),Z)

#endif
