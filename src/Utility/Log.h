/* 
 * File:   Log.h
 * Author: Alberto Lepe <dev@alepe.com>
 *
 * Created on December 1, 2015, 6:00 PM
 */
// Modifications added by JFG

#ifndef LOG_H
#define LOG_H

#include <iostream>

using namespace std;

enum typelog {
    ERROR,  //Uh uh, stuff broke
    WARN,  //Problems arose, but were handled
    INFO,  //Nice to know
    DEBUG, //Ludicrously indepth stuff.
};

struct structlog {
    bool headers = false;
    typelog level = WARN;

    void SetLevel(int i);
};

extern structlog LOGCFG;

class LOG {
public:
    LOG() {}
    LOG(typelog type) {
        msglevel = type;
        if(LOGCFG.headers) {
            operator << (getLabel(type));
        }
    }
    ~LOG() {
        if(opened) {
            clog << "\033[0m" <<endl;
        }
        opened = false;
    }
    template<class T>
    LOG &operator<<(const T &msg) {
        if(msglevel <= LOGCFG.level) {
            clog << msg;
            opened = true;
        }
        return *this;
    }
private:
    bool opened = false;
    typelog msglevel = DEBUG;
    inline string getLabel(typelog type) {
        string label;
        switch(type) {
            case DEBUG: label = "\033[34m[DEBUG] "; break;
            case INFO:  label = "\033[37m[INFO]  "; break;
            case WARN:  label = "\033[33m[WARN]  "; break;
            case ERROR: label = "\033[31m[ERROR] "; break;
        }
        return label;
    }
};

#endif  /* LOG_H */
