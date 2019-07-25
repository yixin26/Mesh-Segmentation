// -----------------------------------------------------------------------------
// ADOBE CONFIDENTIAL
//
// Copyright 2005-2007 Adobe Systems Incorporated
// All Rights Reserved.
// -----------------------------------------------------------------------------
//
// NOTICE:
//    All information  contained  herein is,  and remains the  property of Adobe
// Systems  Incorporated  and  its  supplier,  if  any.   The  intellectual  and
// technical  concepts   contained  herein  are  proprietary  to  Adobe  Systems
// Incorporated  and  its  suppliers  and  may be  covered  by U.S.  and Foreign
// Patents, patents in process,  and are  protected by trade secret or copyright
// law.  Dissemination  of this information  or reproduction of this material is
// strictly forbidden  unless prior  written permission  is obtained  from Adobe
// Systems Incorporated.
// -----------------------------------------------------------------------------


#ifndef __AML_VEC_H
#define __AML_VEC_H

#define USE_GML 0

#if USE_GML

#include "../external/gml_1.1.0/gml.hpp"

namespace AML{
	  typedef double Scalar;
      typedef gml::Vector< 2,Scalar, gml::UNOPT > double2;
      typedef gml::Vector< 3,Scalar, gml::UNOPT > double3;
      typedef gml::Vector< 4,Scalar, gml::UNOPT > double4;
	  typedef gml::Matrix< 2,Scalar, gml::UNOPT > matrix2;
	  typedef gml::Matrix< 3,Scalar, gml::UNOPT > matrix3;
	  typedef gml::Matrix< 4,Scalar, gml::UNOPT > matrix4;

}//namespace AML

#else //#if USE_GML

// Code written by Nathan Carr to represent n-dimensional vector, with 
//   speedups for special case of n=1,2,3. Added to AML by Pushkar Joshi

#include <windows.h>
#include <math.h>

namespace AML {

	typedef double Scalar;


#define GEN_OP_EQUALS(OP,S)\
	inline vec<T,S> &operator OP= (const vec<T,S>& v){\
		for(int i=0;i<S;++i)\
			(*this)[i] OP= v[i];\
		return *this;\
	}\
	inline vec<T,S> &operator OP= (T s){\
		for(int i=0;i<S;++i)\
			(*this)[i] OP= s;\
		return *this;\
	}

#define ALL_VEC_OPS(S) \
	vec<T,S>(){} \
	vec<T,S> (const vec<T,S>& v){\
		for(int i=0;i<S;++i){\
           (*this)[i] = v[i];\
		}\
	}\
	inline vec<T,S> &operator=(const vec<T,S>& v){\
        for(int i=0;i<S;++i)\
            (*this)[i] = v[i];\
        return *this;\
    }\
	inline static unsigned int size(){ return S; } \
	inline void normalize(){T len = length(); if (len) (*this) *= (T)1.0/len; }\
	inline void minOf(const vec<T,S>& v){\
        for(int i=0;i<S;++i)\
			(*this)[i] = ((*this)[i]<v[i])? (*this)[i] : v[i];\
	}\
	inline void maxOf(const vec<T,S>& v){\
        for(int i=0;i<S;++i)\
			(*this)[i] = ((*this)[i]>v[i])? (*this)[i] : v[i];\
	}\
	inline T maxElement(){\
		T result = (*this)[0];\
        for(int i=1;i<S;++i)\
			result = std::max<T>((*this)[i],result);\
		return result;\
	}\
	inline T dot(const vec<T,S>& v)const{\
		T nDot=(*this)[0] * v[0];\
        for(int i=1;i<S;++i)\
			nDot += (*this)[i] * v[i];\
		return nDot;\
	}\
	inline T length()const{\
		return (T)(sqrt( (*this).dot(*this) ));\
	}\
	inline void invert(){\
		for(int i=0;i<S;++i){\
			(*this)[i] = (T)1.0/(*this)[i];\
		}\
	}\
	inline vec<T,S> inverse()const{\
		vec<T,S> result;\
		for(int i=0;i<S;++i){\
			result[i] = 1.0f/(*this)[i];\
		}\
		return result;\
	}\
	inline vec<T,S> operator-()const{\
		vec<T,S> result;\
		for(int i=0;i<S;++i){\
			result[i] = -(*this)[i];\
		}\
		return result;\
	}\
	inline bool operator==(const vec<T,S>& v)const{\
		for(int i=0;i<S;++i){\
			if( (*this)[i] != v[i] )\
				return false;\
		}\
		return true;\
	}\
	inline bool operator!=(const vec<T,S>& v)const{\
		return !( (*this) == v );\
	}\
	GEN_OP_EQUALS(+,S)\
	GEN_OP_EQUALS(-,S)\
	GEN_OP_EQUALS(*,S)\
	GEN_OP_EQUALS(/,S)

#define PARTIAL_VEC_OPS(S) \
	inline T &operator[](size_t nIndex){ return (&x)[nIndex];}\
    inline T operator[](size_t nIndex)const{return (&x)[nIndex];}\
	typedef T* iter;\
	iter begin(){ return &x;}\
	iter end(){   return begin()+S;}

template <class T,unsigned int S> 
class vec
{
public:
    ALL_VEC_OPS(S);
	
	typedef typename T* iter;
	inline iter begin(){ a_data;}
	inline iter end(){   return a_data+S;}
    inline T &operator[](size_t nIndex){ return a_data[nIndex];}
    inline T operator[](size_t nIndex)const{return a_data[nIndex];}
   
private:
    T a_data[S];
};

template <class T> 
class vec<T,2>{
public:
    vec(T s0,T s1):x(s0),y(s1){}
    ALL_VEC_OPS(2);
	PARTIAL_VEC_OPS(2);
	inline void set(T fX,T fY){
		x=fX;y=fY;
	}
	//
    T x,y;
};

template <class T> 
class vec<T,3>{
public:
    vec(T s0,T s1,T s2):x(s0),y(s1),z(s2){}
    ALL_VEC_OPS(3);
	PARTIAL_VEC_OPS(3);
    
	inline void set(T fX,T fY,T fZ){
		x=fX;y=fY;z=fZ;
	}
	inline vec<T,3> cross(const vec<T,3>& v)const{
		return vec<T,3>(y*v.z-z*v.y,
			            z*v.x-x*v.z,
						x*v.y-y*v.x);
	}
	//
    T x,y,z;
};

template <class T> 
class vec<T,4>{
public:
    vec(T s0,T s1,T s2,T s3):x(s0),y(s1),z(s2),w(s3){}
	ALL_VEC_OPS(4);
	PARTIAL_VEC_OPS(4);
    //
    T x,y,z,w;
};


typedef vec<double,2> double2;
typedef vec<double,3> double3;
typedef vec<double,4> double4;
typedef vec<float,2> float2;
typedef vec<float,3> float3;
typedef vec<float,4> float4;
typedef vec<int,2> int2;
typedef vec<int,3> int3;
typedef vec<int,4> int4;
typedef vec<unsigned int,2> uint2;
typedef vec<unsigned int,3> uint3;
typedef vec<unsigned int,4> uint4;
typedef vec<unsigned char,2> ubyte2;
typedef vec<unsigned char,3> ubyte3;
typedef vec<unsigned char,4> ubyte4;
typedef vec<char,2> byte2;
typedef vec<char,3> byte3;
typedef vec<char,4> byte4;

#define BINARY_OP( OP )\
template<class T,unsigned int S> inline \
vec<T,S> operator OP (const vec<T,S>& v0,const vec<T,S>& v1)\
{\
	vec<T,S> result;\
	for(int i=0;i<S;++i)\
       result[i] = v0[i] OP v1[i];\
    return result;\
}\
template<class T,unsigned int S> inline \
vec<T,S> operator OP (const vec<T,S>& v0,T scalar)\
{\
	vec<T,S> result;\
	for(int i=0;i<S;++i)\
       result[i] = v0[i] OP scalar;\
    return result;\
}

BINARY_OP( + )
BINARY_OP( - )
BINARY_OP( * )
BINARY_OP( / )

template<class T,unsigned int S> 
inline vec<T,S> operator * (T scalar,const vec<T,S>& v0)
{
	vec<T,S> result;
	for(int i=0;i<S;++i)
       result[i] = v0[i] * scalar;
    return result;
}

//temp;
template<class T, unsigned int S> inline
bool operator< (const vec<T, S>& v0, const vec<T, S>& v1)
{
	unsigned it = 0;
	while (it<S)
	{
		if (v0[it] != v1[it])
			return v0[it] < v1[it];
		it++;
	}
	return true;
}

} // namespace AML

#endif //#else of #if USE_GML

#endif //__AML_VEC_H