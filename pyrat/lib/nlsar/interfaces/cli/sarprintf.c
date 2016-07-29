#include <stdarg.h>
#include "tools/sarprintf.h"

int sarprintf(const char* format, ...)
{
  va_list args;
  va_start(args, format);
  int res = sarvprintf_std(format, args);
  va_end(args);
  return res;
}

int sarprintf_ret(const char* format, ...)
{
  va_list args;
  va_start(args, format);
  int res = sarvprintf_std_ret(format, args);
  va_end(args);
  return res;
}

int sarprintf_warning(const char* format, ...)
{
  va_list args;
  va_start(args, format);
  int res = sarvprintf_std_warning(format, args);
  va_end(args);
  return res;
}

int sarprintf_error(const char* format, ...)
{
  va_list args;
  va_start(args, format);
  int res = sarvprintf_std_error(format, args);
  va_end(args);
  return res;
}
