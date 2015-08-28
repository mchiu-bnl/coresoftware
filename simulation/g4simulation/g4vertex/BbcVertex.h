#ifndef __BBCVERTEX_H__
#define __BBCVERTEX_H__

#include <phool/PHObject.h>
#include <cmath>
#include <iostream>

class BbcVertex : public PHObject {

public:
  
  virtual ~BbcVertex() {}

  // PHObject virtual overloads
  
  virtual void         identify(std::ostream& os = std::cout) const {os << "BbcVertex base class" << std::endl;}
  virtual BbcVertex*   Clone()                                      {return NULL;}
  virtual void         Reset()                                      {}
  virtual int          IsValid() const                              {return 0;}

  // vertex info
  
  virtual unsigned int get_id() const           {return 0xFFFFFFFF;}
  virtual void         set_id(unsigned int id)  {}
  
  virtual float        get_t0() const           {return NAN;}
  virtual void         set_t0(float t0)         {}

  virtual float        get_t0_err() const       {return NAN;}
  virtual void         set_t0_err(float t0_err) {}
  
  virtual float        get_z() const            {return NAN;}
  virtual void         set_z(float z)           {}

  virtual float        get_z_err() const        {return NAN;}
  virtual void         set_z_err(float z_err)   {}

protected:
  BbcVertex() {}
  
private:
  
  ClassDef(BbcVertex, 1);
};

#endif

