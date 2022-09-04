#ifndef ENUMERATION_H_INCLUDED
#define ENUMERATION_H_INCLUDED

// Enumerates all possible tuples (a_0, a_1, ... a_{n-1})
//    where 0 <= a_i < bounds[i].
struct Enumeration {
    std::vector<int> bounds;
    std::vector<int> v;
    int n;
    bool end;

    Enumeration (int n1)
    {
        n = n1;
        end = false;
        v.resize(n);
        bounds.resize(n);
        for (int i = 0; i < n; ++i)
            bounds[i] = 2;
    }

    Enumeration (std::vector<int> b)
    {
        bounds = b;
        end = false;
        n = b.size();
        v.resize(n);
    }

    void next ()
    {
        bool incremented = false;
        for (int i = 0; i < n; ++i) {
            if (v[i] < bounds[i] - 1) {
                v[i]++;
                for (int j = 0; j < i; ++j)
                    v[j] = 0;
                incremented = true;
                break;
            }
        }
        if (!incremented)
            // for all i=0..n-1: v[i] == bounds[i] - 1
            end = true;
    }
};

#endif // ENUMERATION_H_INCLUDED
