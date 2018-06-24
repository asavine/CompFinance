#pragma once

#include "mcBase.h"
#include "mcMdl.h"
#include "mcPrd.h"
#include <unordered_map>
#include <memory>
using namespace std;

using ModelStore =
unordered_map<int, pair<unique_ptr<Model<double>>, unique_ptr<Model<Number>>>>;
using ProductStore =
unordered_map<int, pair<unique_ptr<Product<double>>, unique_ptr<Product<Number>>>>;

ModelStore modelStore;
ProductStore productStore;

void putBlackScholes(
    const double            spot,
    const double            vol,
    const bool              normal,
    const double            rate,
    const double            div,
    const int               store)
{
    //  We create 2 models, one for valuation and one for risk
    unique_ptr<Model<double>> mdl = make_unique<BlackScholes<double>>(
        spot, vol, normal, rate, div);
    unique_ptr<Model<Number>> riskMdl = make_unique<BlackScholes<Number>>(
        spot, vol, normal, rate, div);

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
    const int               store)
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
const Model<T>* getModel(const int store);

template<>
const Model<double>* getModel(const int store)
{
    auto it = modelStore.find(store);
    if (it == modelStore.end()) return nullptr;
    else return it->second.first.get();
}

template<>
const Model<Number>* getModel(const int store)
{
    auto it = modelStore.find(store);
    if (it == modelStore.end()) return nullptr;
    else return it->second.second.get();
}

void putEuropean(
    const double            strike,
    const Time              exerciseDate,
    const Time              settlementDate,
    const int               store)
{
    //  We create 2 products, one for valuation and one for risk
    unique_ptr<Product<double>> mdl = make_unique<European<double>>(
        strike, exerciseDate, settlementDate);
    unique_ptr<Product<Number>> riskMdl = make_unique<European<Number>>(
        strike, exerciseDate, settlementDate);

    //  And move them into the map
    productStore[store] = make_pair(move(mdl), move(riskMdl));
}

void putBarrier(
    const double            strike,
    const double            barrier,
    const Time              maturity,
    const double            monitorFreq,
    const int               store)
{
    //  We create 2 products, one for valuation and one for risk
    unique_ptr<Product<double>> mdl = make_unique<UOC<double>>(
        strike, barrier, maturity, monitorFreq);
    unique_ptr<Product<Number>> riskMdl = make_unique<UOC<Number>>(
        strike, barrier, maturity, monitorFreq);

    //  And move them into the map
    productStore[store] = make_pair(move(mdl), move(riskMdl));
}

void putEuropeans(
    //  maturities must be given in increasing order
    const vector<double>&   maturities,
    const vector<double>&   strikes,
    const int               store)
{
    //  Create map
    map<Time, vector<double>> options;
    for (size_t i = 0; i < maturities.size(); ++i)
    {
        options[maturities[i]].push_back(strikes[i]);
    }

    //  We create 2 products, one for valuation and one for risk
    unique_ptr<Product<double>> mdl = make_unique<Europeans<double>>(
        options);
    unique_ptr<Product<Number>> riskMdl = make_unique<Europeans<Number>>(
        options);

    //  And move them into the map
    productStore[store] = make_pair(move(mdl), move(riskMdl));
}

template<class T>
const Product<T>* getProduct(const int store);

template<>
const Product<double>* getProduct(const int store)
{
    auto it = productStore.find(store);
    if (it == productStore.end()) return nullptr;
    else return it->second.first.get();
}

template<>
const Product<Number>* getProduct(const int store)
{
    auto it = productStore.find(store);
    if (it == productStore.end()) return nullptr;
    else return it->second.second.get();
}
