#pragma once

#include <pbcopper/utility/Alarm.h>

#define AbortException(msg) PB_ALARM_EXCEPTION_IMPL("pbmm2", msg, "FATAL", "", "AbortException")
