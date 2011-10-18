#include<iostream>
#include<math.h>
#include<fstream>

using namespace std;

//#define EPSILON 0.038776
#define EPSILON 0.019355

//#define NUMBER_STEPS 6
#define NUMBER_STEPS 12

#define NUMBER_SUBSTEPS 2000000


template < typename FunctionType >
double integrate( unsigned int N, double leftB, double rightB, FunctionType& function_to_integrate)
{
    
    if ( rightB < leftB )
    { 
        std :: cout << "Linke Intervallgrenze " << leftB << " kleiner als rechte Intervallgrenze " << rightB << ". Fehler!" << std :: endl;
        abort();
    }

    if ( rightB == leftB )
    { return 0.0;}

    double h = (rightB - leftB) / N;
    
    double integral = 0.0;

    if ( function_to_integrate.antiderivative_available == true )
     { integral = function_to_integrate.evaluate_antiderivative(rightB) - function_to_integrate.evaluate_antiderivative(leftB); }
    else
     {
       for (int i = 0; i < N; i+=1 )
         {
           double x_i = leftB + (i*h);
           integral += h*function_to_integrate.evaluate(x_i);
         }
     }
    return integral;
}





// righ hand side function
class f
{


    public:
    
    double evaluate( const double x )
    {
        return 1.0;
    }

    double evaluate_antiderivative( const double x )
    {
        return x;
    }

    
};




class A
{

    public:

    bool antiderivative_available;

    public:

    A ()
     { antiderivative_available = false; }
    
    double evaluate( const double y )
    {
        //double a_in_y = 2.0 + sin( 2.0 * M_PI * y ); 
        double a_in_y = 1.0 / ( 2.0 + cos( 2.0 * M_PI * y ) );
        return a_in_y;
    }

    double evaluate_antiderivative_of_inverse( const double y )
    {
        double val = 2.0 * y;
        val += ( 1.0 / (2.0 * M_PI ) ) * sin( 2.0 * M_PI * y );
        return val;
    }

    
};


class A_epsilon
{
    public:

    bool antiderivative_available;

    public:
    

    A_epsilon ()
     { antiderivative_available = false; }

    //evaluate a^{\epsilon}
    double evaluate( const double x )
    {
        A diffusion_coefficient;
        const double x_new = x / EPSILON;
        
        double a_eps_in_x = diffusion_coefficient.evaluate( x_new );
        return a_eps_in_x;
        
    }


    //evaluate a^{\epsilon}
    double evaluate_antiderivative_of_inverse( const double x )
    {
        A diffusion_coefficient;
        const double x_new = x / EPSILON;
        double ret = diffusion_coefficient.evaluate_antiderivative_of_inverse( x_new );
        ret *= EPSILON;
        return ret;
    }
    
};


class Exact_Solution
{
};


class Homogenized_Solution
{
};


// (a^{epsilon})^{-1}
class A_epsilon_inverse
{

    public:

    bool antiderivative_available;

public:
    

    A_epsilon_inverse ()
     { antiderivative_available = true; }

    double evaluate( const double x )
    {
        A_epsilon diffusion_coefficient;
        
        double a_eps_in_x = diffusion_coefficient.evaluate( x );
        return (1.0 / a_eps_in_x);
        
    }

    // Stammfunktion (falls bekannt) zum besseren integrieren
    double evaluate_antiderivative( const double x )
    {

//! hard coding!!!!:

        double ret = 2.0*x;
        ret += ( EPSILON / ( 2.0 * M_PI ) ) * sin( 2.0 * ( M_PI / EPSILON ) * x );
        return ret;
        
    }

    
};



// u_n(x) = v_n(x) - g_n(x) ( v_n( x_(n+1) ) / g( x_(n+1) ) )

// v_n( x ) = a^{\eps}( x_{n+1} ) \int_{x_n}^x 1 / ( a^{\eps}(y) ) dy + x_n - x
class V 
{
    
private:
    
    double x_n_;
    double h_;
    
public:
    
    
    V ( double x_n, double h )
    { x_n_ = x_n; h_ = h;}
    
    //evaluate a^{\epsilon}
    double evaluate( const double x )
    {

        A_epsilon a_eps;
        A_epsilon_inverse a_eps_inv;
        
        double x_n_and_1 = x_n_ + h_;
        
        double v_x = 0.0;
        v_x += a_eps.evaluate( x_n_and_1 ) * integrate( NUMBER_SUBSTEPS, x_n_ , x , a_eps_inv );
        v_x -= x - x_n_;  
        
        return v_x;
        
    }
    
    //if g_x is already computed
    inline double evaluate( const double g_x, const double x )
    {
        
        double v_x = 0.0;
        v_x += g_x;
        v_x -= x - x_n_;  
        
        return v_x;
        
    }

};



// u_n(x) = v_n(x) - g_n(x) ( v_n( x_(n+1) ) / g( x_(n+1) ) )

// g_n( x ) = a^{\eps}( x_{n+1} ) \int_{x_n}^x 1 / ( a^{\eps}(y) ) dy
class G 
{
    
private:
    
    double x_n_;
    double h_;
    
public:
    
    
    G ( double x_n, double h )
    { x_n_ = x_n; h_ = h;}
    
    //evaluate a^{\epsilon}
    double evaluate( const double x )
    {
        
        A_epsilon a_eps;
        A_epsilon_inverse a_eps_inv;
        
        double x_n_and_1 = x_n_ + h_;
        
        double g_x = 0.0;
        g_x += a_eps.evaluate( x_n_and_1 ) * integrate( NUMBER_SUBSTEPS, x_n_ , x , a_eps_inv );
        
        return g_x;
        
    }
    
};


int main(int argc, char* argv[])
{
    
    A a;
    A_epsilon a_eps;
    A_epsilon_inverse a_eps_inv;
   
    double left_border = 0.0;
    double right_border = 1.0;
    
    double a_0 = EPSILON / integrate( NUMBER_SUBSTEPS, 0.0, EPSILON, a_eps_inv );
    double a_0_msfem = 0.0;
    
    double h = (right_border - left_border) / NUMBER_STEPS;

    std :: cout << "EPSILON = " << EPSILON << std :: endl;
    std :: cout << "H = " << h << std :: endl << std :: endl;
    std :: cout << "NUMBER_SUBSTEPS = " << NUMBER_SUBSTEPS << std :: endl << std :: endl;

    double msfem_a_0 = 0.0;

    for (int n = 0; n < NUMBER_STEPS; n +=1 )
    {
        
        double x_n = left_border + ( n * h );
        double x_n_and_1 = left_border + ( (n+1) * h );
        
        G g(x_n,h);
        V v(x_n,h);

        double int_a_eps_corrector_n = 0.0;

        double previous_corrector_in_x = 0.0;

        double g_x_n_and_1 = g.evaluate(x_n_and_1);
        // double v_x_n_and_1 = v.evaluate(x_n_and_1);
        double v_x_n_and_1 = v.evaluate(g_x_n_and_1, x_n_and_1);

        for (int sub_n = 0; sub_n < NUMBER_SUBSTEPS; sub_n += 1 )
        {
            // corrector = u_n(x) = v_n(x) - g_n(x) ( v_n( x_(n+1) ) / g( x_(n+1) ) )
            double corrector_in_x = 0.0;
            
            double sub_h = (x_n_and_1 - x_n) / NUMBER_SUBSTEPS;
            double x = x_n + ( sub_h * sub_n );

            double g_x = g.evaluate(x);

            corrector_in_x += v.evaluate( g_x , x );
            corrector_in_x -= ( g_x * v_x_n_and_1 ) / g_x_n_and_1;
            
            double gradient_corrector_in_x = 0.0;
            if ( sub_n != 0 )
             { gradient_corrector_in_x = (corrector_in_x - previous_corrector_in_x) / sub_h; }

            int_a_eps_corrector_n += sub_h * ( gradient_corrector_in_x + 1.0 ) * a_eps.evaluate(x);

            previous_corrector_in_x = corrector_in_x;
        }
        
        int_a_eps_corrector_n /= h;

        a_0_msfem += ( 1.0 / int_a_eps_corrector_n );

        std :: cout << "int_a_eps_corrector_n = " << int_a_eps_corrector_n << std :: endl;
        //std :: cout << "a_0_msfem = " << a_0_msfem << std :: endl;
    }

    a_0_msfem = NUMBER_STEPS * (1.0 / a_0_msfem );

    
    
    
#if 0 
    double integral_a_eps = integrate( NUMBER_STEPS, -1.0, 1.0, a_eps );
    std :: cout << "Integral A_eps = " << integral_a_eps << std :: endl;
   
    double integral_a = integrate( NUMBER_STEPS, -1.0, 1.0, a );
    std :: cout << "Integral A = " << integral_a << std :: endl;

    std :: cout << "corrector_in_x = " << corrector_in_x << std :: endl;  

#endif

    std :: cout << "A^0 = " << a_0 << std :: endl;    
    std :: cout << "A^0_MsFEM = " << a_0_msfem << std :: endl;    
    

    
    return 0;

}
