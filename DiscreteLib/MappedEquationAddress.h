
#pragma once

#include <vector>
#include <map>
#include <algorithm>
#include <cassert>
#include <iostream>


namespace DiscreteLib
{

class IMappedAddress
{
public:
	virtual ~IMappedAddress() {};
    virtual bool hasKey(size_t pt_id) const = 0;
    virtual bool hasValue(size_t eqs_id) const = 0;
    virtual void key_range(size_t &i_start, size_t &i_end) const = 0;
    virtual void activate(size_t pt_id, bool b) = 0;
    virtual bool isActive(size_t pt_id) const = 0;
    virtual void set(size_t pt_id, long eqs_id) = 0;
    virtual size_t setAll(size_t address_start, size_t dn_pt=1) = 0;
    virtual size_t size() const = 0;
    virtual long address(size_t key_id) const = 0;
    virtual long key(long address_id) const = 0;
    virtual bool isSequential() const = 0;
};

class SequentiallyMappedAddress : public IMappedAddress
{
public:
    SequentiallyMappedAddress(size_t pt_id_start, size_t n) : _pt_id_start(pt_id_start), _n(n)
    {
        _dof_begin = 0;
        _delta_per_pt = 1;
    }
    virtual ~SequentiallyMappedAddress() {};

    bool isSequential() const {return true;};

    bool hasKey(size_t pt_id) const
    {
        return (_pt_id_start<=pt_id && pt_id<_pt_id_start+_n);
    }

    bool hasValue(size_t eqs_id) const
    {
        bool in_range = (address(_pt_id_start) <= eqs_id && eqs_id <= address(_pt_id_start+_n-1));
        if (!in_range) return false;
        if ((eqs_id - _dof_begin)%_delta_per_pt!=0) return false;
        return true;
    }

    void key_range(size_t &i_start, size_t &i_end) const
    {
        i_start = _pt_id_start;
        i_end = i_start + _n;
    }

    void activate(size_t pt_id, bool b)
    {
        if (b) {
            _deactive.erase(pt_id);
        } else {
            _deactive.insert(pt_id);
        }
    }
    
    bool isActive(size_t pt_id) const { return _deactive.count(pt_id)==0;};

    void set(size_t pt_id, long eqs_id)
    {
        //invalid
        std::cout << "***Error: called invalid functioons. SequentiallyMappedAddress::set()." << std::endl;
    }

    size_t setAll(size_t address_start, size_t dn_pt)
    {
        _dof_begin = address_start;
        _delta_per_pt = dn_pt;
        return _dof_begin + (_n-_deactive.size())*_delta_per_pt;
    }

    size_t size() const {return _n;};

    long address(size_t pt_id) const
    {
        assert(_pt_id_start<=pt_id && pt_id<_pt_id_start+_n);

        if (_deactive.count(pt_id)>0) return -1;

        if (_pt_id_start<=pt_id && pt_id<_pt_id_start+_n) {
            size_t loc = 0;
            if (_deactive.size()==0) {
                loc = _dof_begin + (pt_id-_pt_id_start)*_delta_per_pt;
            } else {
                size_t inactive_cnt = 0;
                for (size_t i=_pt_id_start; i<pt_id; i++) {
                    if (_deactive.count(i)>0) inactive_cnt++;
                }
                loc = _dof_begin + (pt_id-_pt_id_start-inactive_cnt)*_delta_per_pt;
            }
            return loc;
        } else {
            return -1;
        }
    }

    long key(long address_id) const
    {
        //TODO
        return -1;
    }

private:
    size_t _pt_id_start;
    size_t _n;
    size_t _dof_begin;
    size_t _delta_per_pt;
    std::set<size_t> _deactive;
};

class RandomlyMappedAddress : public IMappedAddress
{
public:
    RandomlyMappedAddress(const std::vector<size_t> &sorted_list_pt_id) : _list_pt_id(sorted_list_pt_id)
    {
        for (size_t i=0; i<sorted_list_pt_id.size(); i++)
            set(sorted_list_pt_id[i], 0);
    }
    RandomlyMappedAddress(size_t pt_id_start, size_t n)
    {
        _list_pt_id.resize(n);
        for (size_t i=0; i<n; i++) {
            _list_pt_id[i] = pt_id_start + i;
            set(_list_pt_id[i], 0);
        }
    }
    virtual ~RandomlyMappedAddress() {};

    bool isSequential() const {return false;};

    bool hasKey(size_t pt_id) const
    {
        return _list_pt_id.end() !=  std::find(_list_pt_id.begin(), _list_pt_id.end(), pt_id);
    }

    bool hasValue(size_t eqs_id) const
    {
        return (_map_pt2eqs.countInB(eqs_id)>0);
    }


    void key_range(size_t &i_start, size_t &i_end) const
    {
        i_start = _list_pt_id[0];
        i_end = _list_pt_id.back()+1;
    }

    void activate(size_t pt_id, bool b)
    {
        if (b) {
            _deactive.erase(pt_id);
        } else {
            _deactive.insert(pt_id);
            set(pt_id, -1);
        }
    }

    bool isActive(size_t pt_id) const { return _deactive.count(pt_id)==0;};

    void set(size_t pt_id, long eqs_id)
    {
        _map_pt2eqs.insert(pt_id, eqs_id);
    }

    size_t setAll(size_t dof_start, size_t delta_per_pt)
    {
        size_t last = 0;
        if (_deactive.size()==0) {
            for (size_t i=0; i<_list_pt_id.size(); i++) {
                last = dof_start + i*delta_per_pt;
                set(_list_pt_id[i], last);
            }
        } else {
            size_t counter = 0;
            for (size_t i=0; i<_list_pt_id.size(); i++) {
                const size_t pt_id = _list_pt_id[i];
                if (_deactive.count(pt_id)>0) continue;
                last = dof_start + counter*delta_per_pt;
                set(pt_id, last);
                counter++;
            }
        }
        return last+1;
    }

    size_t size() const {return _map_pt2eqs.size();};

    long address(size_t pt_id) const
    {
        if (_map_pt2eqs.countInA(pt_id)==0) return -1;
        return _map_pt2eqs.mapAtoB(pt_id);
    }

    long key(long address_id) const
    {
        if (_map_pt2eqs.countInB(address_id)==0) return -1;
        return _map_pt2eqs.mapBtoA(address_id);
    }

private:
    Base::BidirectionalMap<long, long> _map_pt2eqs;
    std::vector<size_t> _list_pt_id;
    std::set<size_t> _deactive;
};

} //end

