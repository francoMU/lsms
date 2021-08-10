
#ifndef MUST_BASE_FUNCTION_HPP
#define MUST_BASE_FUNCTION_HPP


namespace lsms {

    template<class T>
    class BaseFunction {

    public:

        virtual void evaluate(const std::vector <T> &coordinates_vector,
                         T & result,
                         std::vector <T> &gradient) = 0;


    };

}


#endif //MUST_BASE_FUNCTION_HPP
