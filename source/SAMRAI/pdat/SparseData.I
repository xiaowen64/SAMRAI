/********************************************************************** 
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997 - 2011 Lawrence Livermore National Security, LLC
 * Description:   pdat
 **********************************************************************/
namespace SAMRAI {
namespace pdat {
/********************************************************************** 
 * Constructor
 **********************************************************************/ 
template <typename BOX_GEOMETRY>
SAMRAI_INLINE_KEYWORD
SparseData<BOX_GEOMETRY>::Attributes::Attributes(
   const int dsize, const int isize) :
      d_dbl_attrs(dsize), d_int_attrs(isize)
{
}

/********************************************************************** 
 * Copy ctor
 **********************************************************************/ 
template <typename BOX_GEOMETRY>
SAMRAI_INLINE_KEYWORD
SparseData<BOX_GEOMETRY>::Attributes::Attributes(
   const Attributes& other) : 
      d_dbl_attrs(other.d_dbl_attrs),
      d_int_attrs(other.d_int_attrs)
{
}

/********************************************************************** 
 * Dtor
 **********************************************************************/ 
template <typename BOX_GEOMETRY>
SAMRAI_INLINE_KEYWORD
SparseData<BOX_GEOMETRY>::Attributes::~Attributes()
{
}


/********************************************************************** 
 * access operators
 *********************************************************************/  
template <typename BOX_GEOMETRY>
SAMRAI_INLINE_KEYWORD
double&
SparseData<BOX_GEOMETRY>::Attributes::operator[](const DoubleAttributeId& id) 
{
   return d_dbl_attrs[id()];
}

/********************************************************************** 
 *
 *********************************************************************/  
template <typename BOX_GEOMETRY>
SAMRAI_INLINE_KEYWORD
const double&
SparseData<BOX_GEOMETRY>::Attributes::operator[](const DoubleAttributeId& id) const 
{
   return d_dbl_attrs[id()];
}

/********************************************************************** 
 *
 *********************************************************************/  
template <typename BOX_GEOMETRY>
SAMRAI_INLINE_KEYWORD
int&
SparseData<BOX_GEOMETRY>::Attributes::operator[](const IntegerAttributeId& id) 
{
   return d_int_attrs[id()];
}

/********************************************************************** 
 *
 *********************************************************************/  
template <typename BOX_GEOMETRY>
SAMRAI_INLINE_KEYWORD
const int&
SparseData<BOX_GEOMETRY>::Attributes::operator[](const IntegerAttributeId& id) const 
{
   return d_int_attrs[id()];
}

/********************************************************************** 
 * 
 **********************************************************************/ 
template <typename BOX_GEOMETRY>
SAMRAI_INLINE_KEYWORD
bool
SparseData<BOX_GEOMETRY>::Attributes::operator==(
   const Attributes& rhs) const
{
   return d_dbl_attrs == rhs.d_dbl_attrs && 
          d_int_attrs == rhs.d_int_attrs;
}
/********************************************************************** 
 * d'tor
 *********************************************************************/  
template <typename BOX_GEOMETRY>
SAMRAI_INLINE_KEYWORD
SparseData<BOX_GEOMETRY>::~SparseData()
{
}

/********************************************************************** 
 * canEstimateStreamSizeFromBox()
 *********************************************************************/  
template <typename BOX_GEOMETRY>
SAMRAI_INLINE_KEYWORD
bool
SparseData<BOX_GEOMETRY>::canEstimateStreamSizeFromBox() const
{
   return true;
}

/********************************************************************** 
 * empty()
 *********************************************************************/  
template <typename BOX_GEOMETRY>
SAMRAI_INLINE_KEYWORD
bool
SparseData<BOX_GEOMETRY>::empty()
{
   return d_index_to_attribute_map.empty();
}


/********************************************************************** 
 * clear()
 *********************************************************************/  
template <typename BOX_GEOMETRY>
SAMRAI_INLINE_KEYWORD
void
SparseData<BOX_GEOMETRY>::clear()
{
   d_index_to_attribute_map.clear();
}

/********************************************************************** 
 * remove(index)
 *********************************************************************/  
template <typename BOX_GEOMETRY>
SAMRAI_INLINE_KEYWORD
void
SparseData<BOX_GEOMETRY>::remove(Iterator& iterToRemove) 
{
    d_index_to_attribute_map.erase(iterToRemove.d_iterator++);
}

/********************************************************************** 
 * size()
 *********************************************************************/  
template <typename BOX_GEOMETRY>
SAMRAI_INLINE_KEYWORD
int
SparseData<BOX_GEOMETRY>::size() 
{
   return static_cast<int> (d_index_to_attribute_map.size()); 
}

/********************************************************************** 
 * isValid(double_id)
 *********************************************************************/  
template <typename BOX_GEOMETRY>
SAMRAI_INLINE_KEYWORD
bool 
SparseData<BOX_GEOMETRY>::isValid(
   const DoubleAttributeId& id) const
{
   return (id() >= 0) && (id() < d_dbl_attr_size);
}

/********************************************************************** 
 * isValid(int_id)
 *********************************************************************/  
template <typename BOX_GEOMETRY>
SAMRAI_INLINE_KEYWORD
bool 
SparseData<BOX_GEOMETRY>::isValid(
   const IntegerAttributeId& id) const 
{
   return (id() >= 0) && (id() < d_dbl_attr_size);
}


/********************************************************************** 
 * ctor's
 *********************************************************************/  
template <typename BOX_GEOMETRY>
SAMRAI_INLINE_KEYWORD
SparseDataIterator<BOX_GEOMETRY>::SparseDataIterator() :
   d_data(NULL)
{
}

/********************************************************************** 
 * ctor
 *********************************************************************/  
template <typename BOX_GEOMETRY>
SAMRAI_INLINE_KEYWORD
SparseDataIterator<BOX_GEOMETRY>::SparseDataIterator(
   SparseData<BOX_GEOMETRY>& sparse_data) :
      d_data(&sparse_data)
{
   d_iterator = d_data->d_index_to_attribute_map.begin();
}

/********************************************************************** 
 * ctor
 *********************************************************************/  
template <typename BOX_GEOMETRY>
SAMRAI_INLINE_KEYWORD
SparseDataIterator<BOX_GEOMETRY>::SparseDataIterator(
   SparseData<BOX_GEOMETRY>* sparse_data) :
   d_data(sparse_data)
{
   d_iterator = d_data->d_index_to_attribute_map.begin();
}
/********************************************************************** 
 * operator==
 *********************************************************************/  
template <typename BOX_GEOMETRY>
SAMRAI_INLINE_KEYWORD
bool
SparseDataIterator<BOX_GEOMETRY>::operator==(
   const SparseDataIterator<BOX_GEOMETRY>& rhs) const
{
   return d_iterator == rhs.d_iterator;
}

/********************************************************************** 
 * inequality
 *********************************************************************/  
template <typename BOX_GEOMETRY>
SAMRAI_INLINE_KEYWORD
bool
SparseDataIterator<BOX_GEOMETRY>::operator!= (
   const SparseDataIterator<BOX_GEOMETRY>& rhs) const
{
   return ! this->operator==(rhs);
}

/********************************************************************** 
 * getIndex()
 **********************************************************************/ 
template <typename BOX_GEOMETRY>
const hier::Index& 
SparseDataIterator<BOX_GEOMETRY>::getIndex() const
{
   return d_iterator->first;
}
/********************************************************************** 
 * pre-increment operator
 *********************************************************************/  
template <typename BOX_GEOMETRY>
SAMRAI_INLINE_KEYWORD
void
SparseDataIterator<BOX_GEOMETRY>::operator++()
{
   ++d_iterator;
}

/********************************************************************** 
 * post-increment operator
 *********************************************************************/  
template <typename BOX_GEOMETRY>
SAMRAI_INLINE_KEYWORD
void
SparseDataIterator<BOX_GEOMETRY>::operator++(int) {
   d_iterator++;
}

/********************************************************************** 
 * private c'tor
 *********************************************************************/  
template <typename BOX_GEOMETRY>
SAMRAI_INLINE_KEYWORD
SparseDataIterator<BOX_GEOMETRY>::SparseDataIterator(
   SparseData<BOX_GEOMETRY>& sparse_data, 
   typename SparseData<BOX_GEOMETRY>::IndexMap::iterator iterator) :
      d_data(&sparse_data),
      d_iterator(iterator)
{
}


/********************************************************************** 
 * SparseDataAttributeIterator
 *********************************************************************/  
template <typename BOX_GEOMETRY>
SAMRAI_INLINE_KEYWORD
SparseDataAttributeIterator<BOX_GEOMETRY>::SparseDataAttributeIterator(
   const SparseData<BOX_GEOMETRY>& sparse_data,
   const hier::Index& index) :
      d_list(sparse_data._get(index)),
      d_list_iterator(sparse_data._get(index).begin())
{
}
/********************************************************************** 
 * equality operator
 *********************************************************************/  
template <typename BOX_GEOMETRY>
SAMRAI_INLINE_KEYWORD
bool
SparseDataAttributeIterator<BOX_GEOMETRY>::operator==(
   const SparseDataAttributeIterator<BOX_GEOMETRY>& rhs) const
{
   return d_list_iterator == rhs.d_list_iterator;
}

/********************************************************************** 
 * inequality operator
 *********************************************************************/  
template <typename BOX_GEOMETRY>
SAMRAI_INLINE_KEYWORD
bool
SparseDataAttributeIterator<BOX_GEOMETRY>::operator!=(
   const SparseDataAttributeIterator<BOX_GEOMETRY>& rhs) const
{
   return ! this->operator==(rhs);
}

/********************************************************************** 
 * pre-increment operator
 *********************************************************************/  
template <typename BOX_GEOMETRY>
SAMRAI_INLINE_KEYWORD
void
SparseDataAttributeIterator<BOX_GEOMETRY>::operator++()
{
   ++d_list_iterator;
}

/********************************************************************** 
 * post-increment operator
 *********************************************************************/  
template <typename BOX_GEOMETRY>
SAMRAI_INLINE_KEYWORD
void
SparseDataAttributeIterator<BOX_GEOMETRY>::operator++(int)
{
   d_list_iterator++;
}

/********************************************************************** 
 * access operators
 *********************************************************************/  
template <typename BOX_GEOMETRY>
SAMRAI_INLINE_KEYWORD
double&
SparseDataAttributeIterator<BOX_GEOMETRY>::operator[](
   const DoubleAttributeId& id)
{
   return (*d_list_iterator)[id];
}

/********************************************************************** 
 *
 *********************************************************************/  
template <typename BOX_GEOMETRY>
SAMRAI_INLINE_KEYWORD
const double&
SparseDataAttributeIterator<BOX_GEOMETRY>::operator[](
   const DoubleAttributeId& id) const
{
   return (*d_list_iterator)[id];
}

/********************************************************************** 
 *
 *********************************************************************/  
template <typename BOX_GEOMETRY>
SAMRAI_INLINE_KEYWORD
int&
SparseDataAttributeIterator<BOX_GEOMETRY>::operator[](
   const IntegerAttributeId& id)
{
   return (*d_list_iterator)[id];
}

/********************************************************************** 
 *
 *********************************************************************/  
template <typename BOX_GEOMETRY>
SAMRAI_INLINE_KEYWORD
const int&
SparseDataAttributeIterator<BOX_GEOMETRY>::operator[](
   const IntegerAttributeId& id) const
{
   return (*d_list_iterator)[id];
}



}
}
