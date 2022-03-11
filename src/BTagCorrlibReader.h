#ifndef BTagCorrlibReader_H
#define BTagCorrlibReader_H

#include "correction.h"

#include <string>
#include <iostream>

class BTagCorrlibReader
{
    private:
        std::unique_ptr<correction::CorrectionSet> correctionSet_;
    public:
        BTagCorrlibReader()
        {
        }
        
        void loadCorrections(const char* filePath)
        {
            correctionSet_ = std::move(correction::CorrectionSet::from_file(std::string(filePath)));
        }
        
        void printTypeInfo(std::string type)
        {
            if (not correctionSet_) throw std::runtime_error("No btag corrections have been loaded!");
            auto corr = correctionSet_->at(type);
            std::cout<<"info for syst '"<<type<<"'"<<std::endl;
            std::cout<<"  name: "<<corr->name()<<std::endl;
            std::cout<<"  description: "<<corr->description()<<std::endl;
            std::cout<<"  inputs: "<<corr->inputs().size()<<std::endl;
            for (const auto input: corr->inputs())
            {
                std::cout<<"     "<<input.name()<<" ("<<input.typeStr()<<")"<<std::endl;
            }
        }
        
        void printTypes()
        {
            if (not correctionSet_) throw std::runtime_error("No btag corrections have been loaded!");
            std::cout<<"loaded correction types: ";
            for (auto corr: *correctionSet_)
            {
                std::cout<<corr.first<<", ";
            }
            std::cout<<std::endl;
        }
        
        double evaluateBTagWorkingpoint(std::string type, std::string syst, std::string workingpoint, int flavor, float abseta, float pt)
        {
            if (not correctionSet_) throw std::runtime_error("No btag corrections have been loaded!");
            return correctionSet_->at(type)->evaluate({syst,workingpoint,flavor,abseta,pt});
        }

        double evaluateBTagShape(std::string type, std::string syst, int flavor, float abseta, float pt, float discriminant)
        {
            if (not correctionSet_) throw std::runtime_error("No btag corrections have been loaded!");
            return correctionSet_->at(type)->evaluate({syst,flavor,abseta,pt,discriminant});
        }
    
};

#endif
