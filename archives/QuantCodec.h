/** @file QuantCodec.h
 *  @brief This file contains the definitions of the QuantEncoder,
 *         QuantDecoder, QuantEncoderUniform classes.
 *  @author Franck Awounang (awounang)
 *  @bug No known bugs
 */



#ifndef QUANTCODEC_H
#define QUANTCODEC_H
#include <vector>
#include <string>


class QuantEncoder{
public:
    QuantEncoder();
    int getIndex(int v);
    int getValue(int v);
    const std::vector<int>& getLevels();
    ~QuantEncoder(void);
protected:
    std::vector<float> intervals;
    std::vector<int> levels;
    int qMin;
    int qMax;
};


class QuantEncoderUniform: public QuantEncoder {
public:
    QuantEncoderUniform(const std::string& qual, int nbLevels);
    QuantEncoderUniform(const std::vector<int>& qual, int nbLevels);
    QuantEncoderUniform(int qMin, int qMax, int nbLevels);
protected:
    void update(int qMin, int qMax, int nbLevels);
};

class QuantEncoderSemiUniform: public QuantEncoder {
public:
    QuantEncoderSemiUniform(const std::string& qual, int nbLevels);
    QuantEncoderSemiUniform(const std::vector<int>& qual, int nbLevels);
    void update(std::vector<float>& values, std::vector<int>& occurencies, int& dMin, int& dMax, int nbLevels);
};

struct Item {
    int value;
    unsigned int cluster;
    Item(int v, unsigned int c): value(v), cluster(c){}
};

class QuantEncoderLloydMax: public QuantEncoder {
public:
    QuantEncoderLloydMax(const std::string& qual, int nbLevels);
    QuantEncoderLloydMax(const std::vector<int>& data, int nbLevels);
protected:
    void update(std::vector<Item>& items, int nbLevels, int offset, int alphabetSize);
};


class QuantEncoderScore: public QuantEncoder {
public:
    QuantEncoderScore(const std::string& qual, int nbLevels);
    QuantEncoderScore(const std::vector<int>& data, int nbLevels);
protected:
    void update(std::vector<float> pmf, int nbLevels, int offset);
};



class QuantDecoder{
public:
    QuantDecoder(const std::vector<int>& levels);
    int decode(int v);
    ~QuantDecoder(void);
protected:
    std::vector<int> levels;
};



#endif
