#include <string>
#include <vector>
#include <iostream>
#include <sstream>

using std::string;
using std::vector;

struct Options
{
    double p = 1;
    double precision = 6;
    int N = 100;
    int iterationRate = 1000;
    int numOfExperiments = 100;
    bool filtered = false;
};

enum OptionType
{
    OT_P,
    OT_PREC,
    OT_N,
    OT_ITR,
    OT_NOEXP,
    OT_FILT,

    OT_INVALID = -1,
};

vector<string> split(const string text, const char delimiter = '=');
Options getOptions(int argc, char *argv[]);
OptionType getOptionType(string typestr);
bool parseOption(Options *opt, OptionType type, string data);

Options GetOptions(int argc, char *argv[])
{
    Options option;
    OptionType type;
    string message;
    vector<string> splitted;
    int i;

    for (i = 1; i < argc; ++i)
    {
        message = string(argv[i]);
        splitted = split(message);
        type = getOptionType(splitted[0]);
        if (!parseOption(&option, type, splitted[1]))
            std::cout << i << "-th argument is invalid: " << message << std::endl;
    }

    return option;
}

OptionType getOptionType(string typestr)
{
    if (typestr == "p" || typestr == "P")
        return OT_P;
    else if (typestr == "prec" || typestr == "PREC")
        return OT_PREC;
    else if (typestr == "n" || typestr == "N")
        return OT_N;
    else if (typestr == "itr" || typestr == "ITR")
        return OT_ITR;
    else if (typestr == "noexp" || typestr == "NOEXP")
        return OT_NOEXP;
    else if (typestr == "filt" || typestr == "FILTER")
        return OT_FILT;
    else
        return OT_INVALID;
}

bool parseOption(Options *opt, OptionType type, string data)
{
    switch (type)
    {
    case OT_P:
        opt->p = std::stod(data);
        break;
    case OT_PREC:
        opt->precision = std::stoi(data);
        break;
    case OT_N:
        opt->N = std::stoi(data);
        break;
    case OT_ITR:
        opt->iterationRate = std::stoi(data);
        break;
    case OT_NOEXP:
        opt->numOfExperiments = std::stoi(data);
        break;
    case OT_FILT:
        if (data == "true" || data == "TRUE" || data == "t" || data == "T")
            opt->filtered = true;
        else if (data == "false" || data == "FALSE" || data == "f" || data == "F")
            opt->filtered = false;
        else
            return false;
        break;
    default:
        return false;
    }
    return true;
}

vector<string> split(const string text, const char delimiter)
{
    vector<string> columns;

    if (text.empty())
    {
        return columns;
    }

    std::stringstream stream{text};
    string buff;
    while (std::getline(stream, buff, delimiter))
    {
        columns.push_back(buff);
    }
    return columns;
}
