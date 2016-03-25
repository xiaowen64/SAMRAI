//
// File:        main-pointer.C
// Package:     SAMRAI test
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:    $Revision: 586 $
// Modified:    $Date: 2005-08-23 10:49:46 -0700 (Tue, 23 Aug 2005) $
// Description: Main program for testing SAMRAI smart pointers
//

#include "SAMRAI_config.h"

#include "tbox/SAMRAIManager.h"
#include "tbox/MPI.h"
#include "tbox/PIO.h"
#include "tbox/Pointer.h"

using namespace SAMRAI;

class A : public virtual tbox::DescribedClass 
{
public:
   A() {}
   virtual ~A() {}
   virtual void foo(string& trace) = 0;
   virtual void foo1(string& trace) = 0;
};

class B : public virtual tbox::DescribedClass 
{
public:
   B() {}
   virtual ~B() {}
   virtual void foo(string& trace) = 0;
   virtual void foo1(string& trace) { trace = "B::foo1()"; }
};

class C : public A, public B
{
public:
   C() {}
   virtual ~C() {}
   void foo(string& trace) { trace = "C::foo()"; }
   void foo1(string& trace) { trace = "C::foo1()"; }
};

class HaveA 
{
public:
   HaveA(A* a_ptr, 
         const A* a_constptr,
         tbox::Pointer<A> a_smartptr,
         const tbox::Pointer<A> a_constsmartptr) 
   {
      if (a_ptr != a_constptr) {
          tbox::perr << "FAILED: - Test #4a: a_ptr != a_constptr" << endl;
      }
      if (a_ptr != (A*)a_smartptr) {
          tbox::perr << "FAILED: - Test #4b: a_ptr != (A*)a_smartptr" << endl;
      }
      if (a_ptr != (A*)a_constsmartptr) {
          tbox::perr << "FAILED: - Test #4c: a_ptr != (A*)a_constsmartptr" << endl;
      }
      d_a_ptr = a_ptr;
      d_a_smartptr = a_smartptr;
   }
   virtual ~HaveA() {}
   void callFoo()
   {  
      string trace1 = "test-d_a_ptr";
      d_a_ptr->foo(trace1);
      if (trace1 != "C::foo()") {
         tbox::perr << "FAILED: - Test #6a: calling d_a_ptr->foo()" << endl;
         tbox::perr << "   Expected C::foo(), Received " << trace1 << endl;   
      }
      string trace2 = "test-d_a_smartptr";
      d_a_smartptr->foo(trace2);
      if (trace2 != "C::foo()") {
         tbox::perr << "FAILED: - Test #6b: calling d_a_smartptr->foo()" << endl; 
         tbox::perr << "   Expected C::foo(), Received " << trace2 << endl; 
      }
   } 
   void callFoo1()
   {
      string trace = "test-d_a_smartptr";
      d_a_smartptr->foo1(trace);
      if (trace != "C::foo1()") {
         tbox::perr << "FAILED: - Test #8a: calling d_a_smartptr->foo1()" << endl;
         tbox::perr << "   Expected C::foo1(), Received " << trace << endl;
      }
   }
private:
   A* d_a_ptr;
   tbox::Pointer<A> d_a_smartptr;
};

class HaveB 
{
public:
   HaveB(B* b_ptr, 
         const B* b_constptr,
         tbox::Pointer<B> b_smartptr,
         const tbox::Pointer<B> b_constsmartptr) 
   {
      if (b_ptr != b_constptr) {
          tbox::perr << "FAILED: - Test #5b: b_ptr != b_constptr" << endl;
      }
      if (b_ptr != (B*)b_smartptr) {
          tbox::perr << "FAILED: - Test #5b: b_ptr != (A*)b_smartptr" << endl;
      }
      if (b_ptr != (B*)b_constsmartptr) {
          tbox::perr << "FAILED: - Test #5c: b_ptr != (A*)b_constsmartptr" << endl;
      }
      d_b_ptr = b_ptr;
      d_b_smartptr = b_smartptr;
   }    
   virtual ~HaveB() {}
   void callFoo()  
   {
      string trace1 = "test-d_b_ptr";
      d_b_ptr->foo(trace1);
      if (trace1 != "C::foo()") {
         tbox::perr << "FAILED: - Test #7a: calling d_b_ptr->foo()" << endl;
         tbox::perr << "   Expected C::foo(), Received " << trace1 << endl;
      }
      string trace2 = "test-d_b_smartptr";
      d_b_smartptr->foo(trace2);
      if (trace2 != "C::foo()") {
         tbox::perr << "FAILED: - Test #7b: calling d_b_smartptr->foo()" << endl;
         tbox::perr << "   Expected C::foo(), Received " << trace2 << endl;
      }
   } 
   void callFoo1()
   {
      string trace = "test-d_b_smartptr";
      d_b_smartptr->foo1(trace);
      if (trace != "C::foo1()") {
         tbox::perr << "FAILED: - Test #9a: calling d_b_smartptr->foo1()" << endl;
         tbox::perr << "   Expected C::foo1(), Received " << trace << endl;
      }
   }
private:
   B* d_b_ptr;
   tbox::Pointer<B> d_b_smartptr;  
};


class Derived : public tbox::DescribedClass
{
public:
   Derived() { }
   virtual ~Derived() { }
};

class ReallyDerived : public Derived
{
public:
   ReallyDerived() { }
   virtual ~ReallyDerived() { }
};

class ReallyReallyDerived : public ReallyDerived
{
public:
   ReallyReallyDerived() { }
   virtual ~ReallyReallyDerived() { }
};

#include "tbox/ConstPointer.C"
#include "tbox/Pointer.C"
template class tbox::Pointer<A>;
template class tbox::Pointer<B>;
template class tbox::Pointer<C>;
template class tbox::Pointer<Derived>;
template class tbox::Pointer<ReallyDerived>;
template class tbox::Pointer<ReallyReallyDerived>;

/*
 * Function to test casting in function call.
 */

void test(tbox::Pointer<ReallyDerived> a,
          tbox::Pointer<ReallyDerived> b,
          tbox::Pointer<ReallyDerived> c)
{
   if (!a.isNull()) {
      tbox::perr << "FAILED: - Test #2a: in test(), a is non-null" << endl;
   }
   if (b.isNull()) {
      tbox::perr << "FAILED: - Test #2b: in test(), b is null" << endl;
   }
   if (c.isNull()) {
      tbox::perr << "FAILED: - Test #2c: in test(), c is null" << endl;
   }

   tbox::Pointer<ReallyReallyDerived> d = b;
   if (!d.isNull()) {
      tbox::perr << "FAILED: - Test #2d: in test(), d is non-null" << endl;
   }
   tbox::Pointer<ReallyReallyDerived> e = c;
   if (e.isNull()) {
      tbox::perr << "FAILED: - Test #2e: in test(), e is null" << endl;
   }
}

int main( int argc, char *argv[] )
{
   tbox::MPI::init(&argc, &argv);
   tbox::SAMRAIManager::startup();
   tbox::PIO::logAllNodes("pointertest.log");
// tbox::PIO::logOnlyNodeZero("pointertest.log");

   /*
    * Regular pointer tests.
    */

   tbox::Pointer<Derived>       derived        = new Derived;
   tbox::Pointer<ReallyDerived> really_derived = new ReallyDerived;
   tbox::Pointer<ReallyReallyDerived> really_really_derived = 
                                                new ReallyReallyDerived;

   /*
    * Test casting to base class Derived.
    */
   tbox::Pointer<Derived> a = derived;
   if (a.isNull()) {
      tbox::perr << "FAILED: - Test #1a: a is null" << endl;
   }
   tbox::Pointer<Derived> b = really_derived;
   if (b.isNull()) {
      tbox::perr << "FAILED: - Test #1b: b is null" << endl;
   }
   tbox::Pointer<Derived> c = really_really_derived;
   if (c.isNull()) {
      tbox::perr << "FAILED: - Test #1c: c is null" << endl;
   }

   /*
    * Test casting to intermediate class ReallyDerived.
    */
   tbox::Pointer<ReallyDerived> d = derived;
   if (!d.isNull()) {
      tbox::perr << "FAILED: - Test #1d: d is non-null" << endl;
   }
   tbox::Pointer<ReallyDerived> e = really_derived;
   if (e.isNull()) {
      tbox::perr << "FAILED: - Test #1e: e is null" << endl;
   }
   tbox::Pointer<ReallyDerived> f = really_really_derived;
   if (f.isNull()) {
      tbox::perr << "FAILED: - Test #1f: f is null" << endl;
   }

   /*
    * Test casting to most derived class ReallyReallyDerived.
    */
   tbox::Pointer<ReallyReallyDerived> g = derived;
   if (!g.isNull()) {
      tbox::perr << "FAILED: - Test #1g: g is non-null" << endl;
   }
   tbox::Pointer<ReallyReallyDerived> h = really_derived;
   if (!h.isNull()) {
      tbox::perr << "FAILED: - Test #1h: h is non-null" << endl;
   }
   tbox::Pointer<ReallyReallyDerived> i = really_really_derived;
   if (i.isNull()) {
      tbox::perr << "FAILED: - Test #1i: i is null" << endl;
   }

   /*
    * Test casting in function call (#2 tests).
    */
   test(derived, really_derived, really_really_derived);

   
   /*
    * Const pointer tests.
    */

   /*
    * Test casting to base class Derived.
    */ 
   const tbox::Pointer<Derived> j = derived;
   if (j.isNull()) {
      tbox::perr << "FAILED: - Test #3j: j is null" << endl;
   }
   const tbox::Pointer<Derived> k = really_derived;
   if (k.isNull()) {
      tbox::perr << "FAILED: - Test #3k: k is null" << endl;
   }
   const tbox::Pointer<Derived> l = really_really_derived;
   if (l.isNull()) {
      tbox::perr << "FAILED: - Test #3l: l is null" << endl;
   }

   /*
    * Test casting to intermediate class ReallyDerived.
    */
   const tbox::Pointer<ReallyDerived> m = derived;
   if (!m.isNull()) {
      tbox::perr << "FAILED: - Test #3m: m is non-null" << endl;
   }
   const tbox::Pointer<ReallyDerived> n = really_derived;
   if (n.isNull()) {
      tbox::perr << "FAILED: - Test #3n: n is null" << endl;
   }
   const tbox::Pointer<ReallyDerived> o = really_really_derived;
   if (o.isNull()) {
      tbox::perr << "FAILED: - Test #3o: o is null" << endl;
   }

   /*
    * Test casting to most derived class ReallyDerived.
    */
   const tbox::Pointer<ReallyReallyDerived> p = derived;
   if (!p.isNull()) {
      tbox::perr << "FAILED: - Test #3p: p is non-null" << endl;
   }
   const tbox::Pointer<ReallyReallyDerived> q = really_derived;
   if (!q.isNull()) {
      tbox::perr << "FAILED: - Test #3q: q is non-null" << endl;
   }
   const tbox::Pointer<ReallyReallyDerived> r = really_really_derived;
   if (r.isNull()) {
      tbox::perr << "FAILED: - Test #3r: r is null" << endl;
   }

   /*
    * Test casting const pointer to regular pointers.
    */
#if 0
   const tbox::Pointer<ReallyReallyDerived> s = derived;
   if (!s.isNull()) {
      tbox::perr << "FAILED: - Test #3s: s is non-null" << endl;
   }
   const tbox::Pointer<ReallyReallyDerived> t = really_derived;
   if (!t.isNull()) {
      tbox::perr << "FAILED: - Test #3t: t is non-null" << endl;
   }
   const tbox::Pointer<ReallyReallyDerived> u = really_really_derived;
   if (!u.isNull()) {
      tbox::perr << "FAILED: - Test #3u: u is non-null" << endl;
   }
#endif

   /*
    * Virtual base class pointer tests.
    */

   tbox::Pointer<C> my_C = new C();

   // #4 tests
   HaveA my_HaveA(my_C, my_C, my_C, my_C);

   // #5 tests
   HaveB my_HaveB(my_C, my_C, my_C, my_C);

   // #6 tests
   my_HaveA.callFoo();

   // #7 tests
   my_HaveB.callFoo();

   // #8 tests
   my_HaveA.callFoo1();

   // #9 tests
   my_HaveB.callFoo1();

   tbox::pout << "\nPASSED:  testpointer" << endl;

   tbox::SAMRAIManager::shutdown();
   tbox::MPI::finalize();

   return(0); 
}
