#include <sys/time.h>
//int gettime_(null){
int gettime_(){
	struct timeval tp;
	struct timezone tzp;

	gettimeofday(&tp,&tzp);
	return (tp.tv_usec);
}
