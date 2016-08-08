/** @file debug.h
 *  @brief This file contains useful debug macros.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  YYYY-MM-DD: What (who)
 */

#ifndef DEBUG_H
#define DEBUG_H

#include "Common/fileSystemHelpers.h"

// Safe debug macro, usage example: DEBUG("status: %d\n", s);
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

// This expands to the module name, usage example: std::cout << ME << "Hello";
#define MODULE_NAME removeFileNameExtension(fileBaseName(std::string(__FILE__)))
#define COUT_PREFIX std::string("[" + MODULE_NAME + "] ")
#define ME COUT_PREFIX // for brevity

#endif // DEBUG_H

