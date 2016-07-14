/** @file common.h
 *  @brief This file contains common global definitions.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  2017-07-14: Added functions fileBaseName and removeFileNameExtension (voges)
 *  2016-07-14: Added module name macro (voges)
 *  2016-05-11: Added safe debug macro (voges)
 */

#ifndef COMMON_H
#define COMMON_H

#include <string>

// safe debug macro, usage example: DEBUG("status: %d\n", s);
#define DBG
#ifdef DBG
    #define DEBUG(c,...)\
        do {\
            fprintf(stderr, "%s:%s:%d: " c, __FILE__, __FUNCTION__, \
                    __LINE__, ##__VA_ARGS__); \
        } while (false)
#else
    #define DEBUG(c,...) do { } while (false)
#endif

// this expands to the module name, usage example: std::cout << ME << "Hello";
#define MODULE_NAME removeFileNameExtension(fileBaseName(std::string(__FILE__)))
#define COUT_PREFIX std::string("[" + MODULE_NAME + "] ")
#define ME COUT_PREFIX // for brevity

bool fileExists(const std::string &fileName);
std::string fileBaseName(const std::string &path);
std::string fileNameExtension(const std::string &fileName);
std::string removeFileNameExtension(const std::string &fileName);
bool yesno(void);

#endif // COMMON_H

