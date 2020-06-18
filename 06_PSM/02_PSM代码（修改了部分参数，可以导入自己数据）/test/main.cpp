/***************************************************************************
                         main.cpp  -  Starts the scan matching test.
                           -------------------
   begin                : July 2010
   version              : 0.1
   copyright            : (C) 2010 by Albert Diosi
   email                : albert.diosi@gmail.com
***************************************************************************/
/****************************************************************************
Copyright (c) 2010-2015, Albert Diosi
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
    * The name of the copyright holders may not be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
****************************************************************************/
#include <string.h>

#include "polar_match.h"

/** @brief Starts the scan matching unit test.
*/
int main(int argc, char *argv[])
{
  printf("Polar Scan Matching tests. Run with --non-interactive for no graphics.测试输出1\n");
  bool interactive = true;

  if(argc==2)
  {
    if(!strcmp("--non-interactive", argv[1]))
      interactive = false;
  }
  printf("测试输出2\n");
  pm_unit_test(PM_PSM, interactive);
  //pm_unit_test(PM_ICP, interactive);
  return 0;
}
