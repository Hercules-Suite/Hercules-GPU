#include <mpi.h>

enum TimerKind { UNDEF_TIMERKIND=0, NON_CUMULATIVE=1, MAX=2, MIN=4, AVERAGE=8 };

extern void     Timer_Start(const char* TimerName);
extern void     Timer_Stop(const char* TimerName);
extern void     Timer_Reset(const char* TimerName);
extern int      Timer_Exists(const char* TimerName);
extern void     Timer_Config(const char* TimerName, enum TimerKind kind);
extern double   Timer_Value(const char* TimerName, enum TimerKind kind);
extern void     Timer_PrintAll(MPI_Comm Communicator);
extern void     Timer_Reduce(const char* TimerName, enum TimerKind kind, MPI_Comm Communicator);

