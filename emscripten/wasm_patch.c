#include <pthread.h>

static int dummy_thread_struct = 0;

int pthread_create(pthread_t *thread, const pthread_attr_t *attr, void *(*start_routine) (void *), void *arg) {
    // Trick C++ into thinking the thread started successfully
    if (thread) *thread = (pthread_t)&dummy_thread_struct;
    return 0;
}

int pthread_detach(pthread_t thread) { return 0; }
int pthread_join(pthread_t thread, void **retval) { return 0; }