#ifndef UTIL_HPP
#define UTIL_HPP


void handler(int sig);
void __terminal_error(const char *file, const unsigned line, const char *func, const char *msg);
inline int rand_int(int max_int);
inline int rand_int_range(int min_int, int max_int);
int count_lines(FILE *fp);
void write_ints_to_fp(FILE *fp, int *data, int N);
void write_ints_to_file(const char *fn, int *data, int N);
void write_flt_vec2_to_file(const char *fn, vector<vector<float>> *vec2);
void write_dbl_vec2_to_file(const char *fn, vector<vector<double>> *vec2);
std::istream &safeGetline(std::istream &is, std::string &t); 
int *getRandomOrderInts(int rangeEnd); 

#include <cstdlib> // exit




/*#define STATCODE(a) a*/
#define STATCODE(a)
#define terminal_error(errstr) __terminal_error(__FILE__, __LINE__, __func__, errstr)
#define safemalloc(n) safe_malloc(n, __LINE__)

#endif
