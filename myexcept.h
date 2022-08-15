#ifndef MYEXCEPT_H
#define MYEXCEPT_H
#include<vector>
using namespace std;
namespace MyMath {
    class MyExcept
    {
        
    public:
        MyExcept(string msg)
        {
            m_msg = msg;
        }
        string getMsg()
        {
            return m_msg;
        }
    private:
        string m_msg;
        
    };//MyExcept
}//MyMath
#endif