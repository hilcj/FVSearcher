#ifndef FVS_COMMON_H_
#define FVS_COMMON_H_
// xyz to flat id
#include <string>
using namespace std;

// ijk to flat id
inline void ijk2flat(int il, int jl, int kl, int i, int j, int k, int & flat){
	flat = i * jl * kl + j * kl + k;
}

inline void ijk2flat(int il, int jl, int kl, int i, int j, int k, long long & flat){
    flat = i * jl * kl + j * kl + k;
}

// flat id to ijk
inline void flat2ijk(int il, int jl, int kl, int & i, int & j, int & k, int flat){
	k = flat % kl;
	j = (flat / kl) % jl;
	i = (flat / kl / jl);
}

inline void flat2ijk(int il, int jl, int kl, int & i, int & j, int & k, long long flat){
    k = flat % kl;
    j = (flat / kl) % jl;
    i = (flat / kl / jl);
}

// Trim space at the begin and end of string
string trim_string(const string & str){
    size_t first = str.find_first_not_of(' ');
    if (string::npos == first){
        return str;
    }
    size_t last = str.find_last_not_of(' ');
    return str.substr(first, (last - first + 1));
}

// Check if string ends with given substring. 
bool ends_with(string const & value, string const & ending)
{
    if (ending.size() > value.size()) return false;
    return equal(ending.rbegin(), ending.rend(), value.rbegin());
}
#endif
