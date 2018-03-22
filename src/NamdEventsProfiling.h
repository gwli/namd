/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef NAMDEVENTSPROFILING_H

#define NAMDEVENTSPROFILING_H

#include "common.h"
#include <map>

struct NamdProfileEvent {
  typedef enum {
    #define NAMD_PROFILE_EVENT(a, b) a,
    #include "NamdEventsProfiling.def"
    #undef NAMD_PROFILE_EVENT
    EventsCount
  } Event;
};

char const* const NamdProfileEventStr[] = {
  #define NAMD_PROFILE_EVENT(a,b) b,
  #include "NamdEventsProfiling.def"
  #undef NAMD_PROFILE_EVENT
  0
};

//
// When defined, use NVTX to record CPU activity into CUDA profiling run.
//
#undef NAMD_EVENT_START
#undef NAMD_EVENT_STOP
#undef NAMD_EVENT_RANGE

#if defined(NAMD_CUDA) && defined(NAMD_USE_NVTX)

#include <nvToolsExt.h>

// C++ note: declaring const variables implies static (internal) linkage,
// and you have to explicitly specify "extern" to get external linkage.
const uint32_t NAMD_nvtx_colors[] = {
  0x0000ff00,
  0x000000ff,
  0x00ffff00,
  0x00ff00ff,
  0x0000ffff,
  0x00ff0000,
  0x00ffffff,
};
const int NAMD_nvtx_colors_len = sizeof(NAMD_nvtx_colors)/sizeof(uint32_t);

#define NAMD_REGISTER_EVENT(name,cid) do { } while(0) // must terminate with semi-colon

// start recording an event
#define NAMD_EVENT_START(eon,id) \
  do { \
    if (eon) { \
      int color_id = id; \
      color_id = color_id % NAMD_nvtx_colors_len; \
      nvtxEventAttributes_t eventAttrib = {0}; \
      eventAttrib.version = NVTX_VERSION; \
      eventAttrib.size = NVTX_EVENT_ATTRIB_STRUCT_SIZE; \
      eventAttrib.colorType = NVTX_COLOR_ARGB; \
      eventAttrib.color = NAMD_nvtx_colors[color_id]; \
      eventAttrib.messageType = NVTX_MESSAGE_TYPE_ASCII; \
      eventAttrib.message.ascii = NamdProfileEventStr[id]; \
      nvtxRangePushEx(&eventAttrib); \
    } \
  } while(0)  // must terminate with semi-colon

// stop recording an event
#define NAMD_EVENT_STOP(eon,id) \
  do { \
    if (eon) { \
      nvtxRangePop(); \
    } \
  } while (0)  // must terminate with semi-colon

// embed event recording in class to automatically pop when destroyed
class NAMD_NVTX_Tracer {
  protected:
    int evon;  // is event on?
  public:
    NAMD_NVTX_Tracer(int eon, int id = 0) : evon(eon) {
      NAMD_EVENT_START(eon, id);
    }
    ~NAMD_NVTX_Tracer() { NAMD_EVENT_STOP(evon, 0); }
};

// call NAMD_EVENT_RANGE at beginning of function to push event recording
// destructor is automatically called on return to pop event recording
#define NAMD_EVENT_RANGE(eon,id) \
  NAMD_NVTX_Tracer namd_nvtx_tracer(eon,id)
  // must terminate with semi-colon

#elif CMK_TRACE_ENABLED

//
// Otherwise use Projections user events for profiling if Projections tracing is enabled
//
#define SEQUENCER_EVENT_ID_START 150
//TODO: Remove after support is added on Charm/Projections internals
#define traceUserBracketEventStart(eid)  do { } while(0)
#define traceUserBracketEventStop(eid)   do { } while(0)

#define NAMD_REGISTER_EVENT(name,id) \
  do { \
    int eventID = SEQUENCER_EVENT_ID_START+id; \
    traceRegisterUserEvent(name, eventID); \
  } while(0) // must terminate with semi-colon

#define NAMD_EVENT_START(eon,id) \
  do {\
    if (eon) { \
      int eventID = SEQUENCER_EVENT_ID_START+id; \
      traceBeginUserBracketEvent(eventID); \
    } \
  } while(0) // must terminate with semi-colon

#define NAMD_EVENT_STOP(eon,id) \
  do { \
    if (eon) { \
      int eventID = SEQUENCER_EVENT_ID_START+id; \
      traceEndUserBracketEvent(eventID); \
    } \
  } while(0)  // must terminate with semi-colon

class NAMD_Sequencer_Events_Tracer {
  int tEventID;
  int tEventOn;
  public:
    NAMD_Sequencer_Events_Tracer(int eon, int id = 0)
      : tEventOn(eon), tEventID(id) {
      NAMD_EVENT_START(tEventOn, tEventID);
    }
    ~NAMD_Sequencer_Events_Tracer() {
      NAMD_EVENT_STOP(tEventOn, tEventID);
    }
};

// call NAMD_EVENT_RANGE at beginning of function to push event recording
// destructor is automatically called on return to pop event recording
#define NAMD_EVENT_RANGE(eon,id) \
    NAMD_Sequencer_Events_Tracer namd_events_tracer(eon,id)
    // must terminate with semi-colon

#else

//
// Otherwise all profiling macros become no-ops.
//
#define NAMD_REGISTER_EVENT(name,cid)     do { } while(0)  // must terminate with semi-colon
#define NAMD_EVENT_START(eon,id)          do { } while(0)  // must terminate with semi-colon
#define NAMD_EVENT_STOP(eon,id)           do { } while(0)  // must terminate with semi-colon
#define NAMD_EVENT_RANGE(eon,id)          do { } while(0)  // must terminate with semi-colon

#endif // NAMD_CUDA && NAMD_USE_NVTX

#endif /* NAMDEVENTSPROFILING_H */
