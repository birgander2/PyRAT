#include "tools/sarwaitbar.h"

int sarwaitbar_open(void)
{
  return sarwaitbar_std_open();
}

int sarwaitbar_update(int percent)
{
  return sarwaitbar_std_update(percent);
}

int sarwaitbar_close(void)
{
  return sarwaitbar_std_close();
}
