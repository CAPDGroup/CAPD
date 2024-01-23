/*
 * File:   removeUnused.h
 * Author: kapela
 *
 * Created on 17 luty 2015, 02:55
 */

#ifndef CAPD_AUXIL_IGNOREUNUSED_H
#define	CAPD_AUXIL_IGNOREUNUSED_H

namespace capd{
  namespace auxil{

    template <class T1>
    inline void ignoreUnused(const T1&) {}

    template <class T1, class T2>
    inline void ignoreUnused(const T1&, const T2&) {}

    template <class T1, class T2, class T3>
    inline void ignoreUnused(const T1&, const T2&, const T3&) {}

    template <class T1, class T2, class T3, class T4>
    inline void ignoreUnused(const T1&, const T2&, const T3&, const T4&) {}

    template <class T1, class T2, class T3, class T4, class T5>
    inline void ignoreUnused(const T1&, const T2&, const T3&, const T4&, const T5&) {}


  }
}

#endif	/* CAPD_AUXIL_IGNOREUNUSED_H */
