
/*
Written by Antoine Savine in 2018

This code is the strict IP of Antoine Savine

License to use and alter this code for personal and commercial applications
is freely granted to any person or company who purchased a copy of the book

Modern Computational Finance: AAD and Parallel Simulations
Antoine Savine
Wiley, 2018

As long as this comment is preserved at the top of the file
*/

#pragma once

//  Memory storage for models and products,
//  See chapter 6

#include "mcBase.h"
#include "mcMdl.h"
#include "mcPrd.h"
#include <unordered_map>
#include <memory>
using namespace std;

using ModelStore =
unordered_map<string, pair<unique_ptr<Model<double>>, unique_ptr<Model<Number>>>>;
using ProductStore =
unordered_map<string, pair<unique_ptr<Product<double>>, unique_ptr<Product<Number>>>>;

ModelStore modelStore;
ProductStore productStore;

void putBlackScholes(
    const double            spot,
    const double            vol,
    const bool              qSpot,
    const double            rate,
    const double            div,
    const string&           store)
{
    //  We create 2 models, one for valuation and one for risk
    unique_ptr<Model<double>> mdl = make_unique<BlackScholes<double>>(
        spot, vol, qSpot, rate, div);
    unique_ptr<Model<Number>> riskMdl = make_unique<BlackScholes<Number>>(
        spot, vol, qSpot, rate, div);

    //  And move them into the map
    modelStore[store] = make_pair(move(mdl), move(riskMdl));
}

void putDupire(
    const double            spot,
    const vector<double>&   spots,
    const vector<Time>&     times,
    //  spot major
    const matrix<double>&   vols,
    const double            maxDt,
    const string&           store)
{
    //  We create 2 models, one for valuation and one for risk
    unique_ptr<Model<double>> mdl = make_unique<Dupire<double>>(
        spot, spots, times, vols, maxDt);
    unique_ptr<Model<Number>> riskMdl = make_unique<Dupire<Number>>(
        spot, spots, times, vols, maxDt);

    //  And move them into the map
    modelStore[store] = make_pair(move(mdl), move(riskMdl));
}

template<class T>
const Model<T>* getModel(const string& store);

template<>
const Model<double>* getModel(const string& store)
{
    auto it = modelStore.find(store);
    if (it == modelStore.end()) return nullptr;
    else return it->second.first.get();
}

template<>
const Model<Number>* getModel(const string& store)
{
    auto it = modelStore.find(store);
    if (it == modelStore.end()) return nullptr;
    else return it->second.second.get();
}

pair<const vector<string>*, const vector<double*>*> getModelParameters(const string& store)
{
    auto it = modelStore.find(store);
    if (it == modelStore.end()) return make_pair(nullptr, nullptr);
    else
    {
        auto* mdl = it->second.first.get();
        return make_pair(&mdl->parameterLabels(),&mdl->parameters());
    }
}

void putEuropean(
    const double            strike,
    const Time              exerciseDate,
    const Time              settlementDate,
    const string&           store)
{
    //  We create 2 products, one for valuation and one for risk
    unique_ptr<Product<double>> prd = make_unique<European<double>>(
        strike, exerciseDate, settlementDate);
    unique_ptr<Product<Number>> riskPrd = make_unique<European<Number>>(
        strike, exerciseDate, settlementDate);

    //  And move them into the map
    productStore[store] = make_pair(move(prd), move(riskPrd));
}

void putBarrier(
    const double            strike,
    const double            barrier,
    const Time              maturity,
    const double            monitorFreq,
    const double            smooth,
    const string&           store)
{
    const double smoothFactor = smooth <= 0 ? EPS : smooth;

    //  We create 2 products, one for valuation and one for risk
    unique_ptr<Product<double>> prd = make_unique<UOC<double>>(
        strike, barrier, maturity, monitorFreq, smoothFactor);
    unique_ptr<Product<Number>> riskPrd = make_unique<UOC<Number>>(
        strike, barrier, maturity, monitorFreq, smoothFactor);

    //  And move them into the map
    productStore[store] = make_pair(move(prd), move(riskPrd));
}

void putContingent(
    const double            coupon,
    const Time              maturity,
    const double            payFreq,
    const double            smooth,
    const string&           store)
{
    const double smoothFactor = smooth <= 0 ? 0.0 : smooth;

    //  We create 2 products, one for valuation and one for risk
    unique_ptr<Product<double>> prd = make_unique<ContingentBond<double>>(
        maturity, coupon, payFreq, smoothFactor);
    unique_ptr<Product<Number>> riskPrd = make_unique<ContingentBond<Number>>(
        maturity, coupon, payFreq, smoothFactor);

    //  And move them into the map
    productStore[store] = make_pair(move(prd), move(riskPrd));
}

void putEuropeans(
    //  maturities must be given in increasing order
    const vector<double>&   maturities,
    const vector<double>&   strikes,
    const string&           store)
{
    //  Create map
    map<Time, vector<double>> options;
    for (size_t i = 0; i < maturities.size(); ++i)
    {
        options[maturities[i]].push_back(strikes[i]);
    }

    //  We create 2 products, one for valuation and one for risk
    unique_ptr<Product<double>> prd = make_unique<Europeans<double>>(
        options);
    unique_ptr<Product<Number>> riskPrd = make_unique<Europeans<Number>>(
        options);

    //  And move them into the map
    productStore[store] = make_pair(move(prd), move(riskPrd));
}

template<class T>
const Product<T>* getProduct(const string& store);

template<>
const Product<double>* getProduct(const string& store)
{
    auto it = productStore.find(store);
    if (it == productStore.end()) return nullptr;
    else return it->second.first.get();
}

template<>
const Product<Number>* getProduct(const string& store)
{
    auto it = productStore.find(store);
    if (it == productStore.end()) return nullptr;
    else return it->second.second.get();
}

const vector<string>* getPayoffLabels(const string& store)
{
    auto it = productStore.find(store);
    if (it == productStore.end()) return nullptr;
    else return & it->second.first.get()->payoffLabels();

}