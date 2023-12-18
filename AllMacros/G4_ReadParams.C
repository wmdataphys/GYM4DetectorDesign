#ifndef MACRO_G4READPARAMS_C
#define MACRO_G4READPARAMS_C

#include <string>
#include <map>
#include <TString.h>
#include <cctype>

int ReadParams(string FileName, map <string, double> *Params)
{
    /*
    ReadParams(FileName) reads the parameters from the file FileName.
    The file FileName must contain params in the following order
    (one param per line):
    "param1" : value1
    "param2" : value2
    ...
    */
    std::ifstream ParamsFile(FileName.c_str());
    if (!ParamsFile.is_open())
    {
        std::cout << "Error opening file " << FileName << endl;
        return 1;
    }
    string Line;
    while(std::getline(ParamsFile, Line))
    {
        int isAlpha = isalpha(Line[0]);
        if (isAlpha != 0) {
        string ParamName = Line.substr(0, Line.find(" : "));
        string ParamValue = Line.substr(Line.find(":") + 1, Line.length());
        //std::cout << Line <<  ", ParamValue Name : " << ParamName << ", ParamValue = " << ParamValue << std::endl;
        //ParamName.erase(std::remove_if(ParamName.begin(), ParamName.end(), " "), ParamName.end());
        (*Params)[ParamName] = stod(ParamValue);
        }
    }
    ParamsFile.close();

    /*
    for(auto it = Params->begin(); it != Params->end(); ++it)
    {
        std::cout << it->first << " -- " << it->second << std::endl;
    }
    std::cout << "Total number of Parameters found = " << Params->size() << std::endl;
    */
    return 0;

}

#endif /* MACRO_G4READPARAMS_C */