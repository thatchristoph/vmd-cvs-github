
#ifndef HASHARRAY_H
#define HASHARRAY_H

struct hasharray;
typedef struct hasharray hasharray;

hasharray * hasharray_create(void **itemarray, int itemsize);
void hasharray_destroy(hasharray *a);

int hasharray_insert(hasharray *a, const char *key);

int hasharray_delete(hasharray *a, const char *key);

#define HASHARRAY_FAIL -1

int hasharray_index(hasharray *a, const char *key);

int hasharray_count(hasharray *a);

#endif

