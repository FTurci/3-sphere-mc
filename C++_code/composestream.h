#ifndef __COMPOSESTREAM_H
#define __COMPOSESTREAM_H 

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <functional>

class ComposeStream: public std::ostream
{
    struct ComposeBuffer: public std::streambuf
    {
        void addBuffer(std::streambuf* buf)
        {
            bufs.push_back(buf);
        }
        virtual int overflow(int c)
        {
            std::for_each(bufs.begin(),bufs.end(),std::bind2nd(std::mem_fun(&std::streambuf::sputc),c));
            return c;
        }

        private:
            std::vector<std::streambuf*>    bufs;

    };  
    ComposeBuffer myBuffer;
    public: 
        ComposeStream()
            :std::ostream(NULL)
        {
            std::ostream::rdbuf(&myBuffer);
        }   
        void linkStream(std::ostream& out)
        {
            out.flush();
            myBuffer.addBuffer(out.rdbuf());
        }
};

#endif