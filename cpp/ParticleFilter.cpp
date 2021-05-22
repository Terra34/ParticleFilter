#include <iostream>
#include <random>
#include <math.h>

#define PI 3.14159265359

template<std::size_t N>
class ParticleFilter {
    public:
        float filterState[2];
        float t;
        int numParticles = N;

        ParticleFilter(float pt = 1.0) : filterState{0.0, 0.0} 
        {
            /*
                pt => time parameter
            */
            t = pt;
            setGenerator();
            std::uniform_real_distribution<float> resampleRandom(0, 1.);
            std::uniform_real_distribution<float> distribution(-180, +180);
            for (int i = 0; i < numParticles; i++){
                weights[i] = 1./numParticles;
                for (int j = 0; j < 2; j++){
                    _particles[i][j] = distribution(generator);
                }
            }
            particles = _particles;
            copyParticles = _copyParticles;
        }

        void estimate() 
        {
            /*
                set estimated state
            */
            filterState[0] = 0.;
            filterState[1] = 0.;
            for (int i = 0; i < numParticles; i++){
                for (int j = 0; j < 2; j++){
                    filterState[j] += particles[i][j] * weights[i] / sum_weights;
                }
            }
            return;
        }

        void predict(float* gm, float pstd = 1.0) 
        {
            /*
                gm => gyroscope measurements
                pstd => random noise standard deviation

                predict next state of particles using state model
            */
            std::normal_distribution<float> distro(0.0, pstd);

            for (int i = 0; i < numParticles; i++){
                particles[i][0] += *gm * t + distro(generator);
                particles[i][1] += *(gm+1) * t + distro(generator);

                // check if value got out of [-180, 180] range
                if (particles[i][0] > 180)
                    particles[i][0] -= 360;
                else if (particles[i][0] < -180)
                    particles[i][0] += 360;

                if (particles[i][1] > 180)
                    particles[i][1] -= 360;
                else if (particles[i][1] < -180)
                    particles[i][1] += 360;
            }
            return;
        }

        void update(float* mS, float* sE)
        {
            /*  
                mS => measured state
                sE => sensor errors

                update weights based on similarity
            */
            sum_weights = 0.0;
            for (int i=0; i<numParticles; i++){
                weights[i] *= 1. / ( pow(particles[i][0] - *mS, 2) + pow(particles[i][1] - *(mS + 1), 2) );
                /*
                for (int j=0; j<2; j++){
                    weights[i] *= normProb(particles[i][j], *(mS + j), *(sE + j));
                }
                */
                // add float smallest value (no division by 0 errors)
                weights[i] += std::numeric_limits<float>::min();
                sum_weights += weights[i];
            }
            return;
        }

        void resample()
        {
            /*
                resample particles based on weights
            */
            int indexes[numParticles];
            float residual;
            float residuals[numParticles];
            float sum_residual=0;
            int k = 0;
            int num_copies;
            float w;
            float random;

            // initial sampling and calculating residuals
            for (int i=0; i<numParticles; i++){
                w = numParticles * weights[i] / sum_weights;
                num_copies = (int)w;
                residual = w - num_copies;
                residuals[i] = residual;
                sum_residual += residual;

                if (num_copies > 0)
                    for (int j=0; j<num_copies; j++)
                        indexes[k++] = i;
            }
            
            // cumulative sum of residuals
            residuals[0] = residuals[0] / sum_residual;
            for (int i=1; i<numParticles; i++){
                residuals[i] = residuals[i] / sum_residual;
                residuals[i] += residuals[i-1];
            }
            // ensure it's a real CDF
            residuals[numParticles-1] = 1.;

            // sample from residuals
            while (k < numParticles){
                random = resampleRandom(generator);
                indexes[k++] = searchSorted(residuals, random);
            }

            resampleParticles(indexes);
            return;
        }

    private:
        float _particles[N][2];
        float _copyParticles[N][2];
        float (*particles)[2], (*copyParticles)[2];
        float weights[N];
        float sum_weights = 1.;
        std::default_random_engine generator;
        std::uniform_real_distribution<float> resampleRandom;

        void setGenerator()
        {
            std::random_device myRandomDevice;
            unsigned seed = myRandomDevice();
            std::default_random_engine gen(seed);
            generator = gen;
            return;
        }

        int searchSorted(float* arr, float r)
        {
            /*
                arr => Array that needs to be searched
                r => Value to search by
            */
            int L = 0;
            int R = numParticles - 1;
            int m;

            while (R >= L)
            {
                m = (int) ((L+R) / 2);
                if (r > *(arr + m))
                {
                    if (r <= *(arr + m + 1) || (m + 1 > numParticles - 1))
                        return m+1;
                    L = m + 1;
                }
                else // if (r <= *(arr + m))
                {
                    if (r > *(arr + m - 1) || (m - 1 < 0))
                        return m;
                    R = m - 1;
                }
            }
            return -1;
        }

        float normProb(float x, float mean, float std)
        {
            /*
                p(x|mu, std)
            */
            return (1 / std / sqrt(2*PI)) * exp(-1 * pow(x - mean, 2) / (2 * pow(std, 2)));
        }  

        void resampleParticles(int * indexes)
        {
            float (*temp)[2];
            for(int i=0; i<numParticles; i++){
                for(int j=0; j<2; j++)
                    copyParticles[i][j] = particles[*(indexes + i)][j];
                weights[i] = 1./numParticles;
            }
            temp = particles;
            particles = copyParticles;
            copyParticles = temp;
            return;
        } 
};

float * getPitchRoll(float* meas)
{
    // meas is acc_data => [accx accy accz]
    static float pitchRoll[2];
    // roll
    pitchRoll[0] = atan2(-1. * *(meas + 1), sqrt(pow(*meas, 2) + pow(*(meas + 2), 2))) * 180 / PI;

    // pitch
    pitchRoll[1] = atan2(*meas, sqrt(pow(*(meas + 1), 2) + pow(*(meas + 2), 2))) * 180 / PI;

    return pitchRoll;
}

float getYaw(float* fS, float* meas)
{
    /*  
        fS => filtered state - assumes degrees!
        meas => magnetometer data [magx magy magz]
    */
   float sina, sinb, cosa, cosb, XH, YH;

   sina = sin(*fS * PI / 180);
   sinb = sin(*(fS + 1) * PI / 180);
   cosa = cos(*fS * PI / 180);
   cosb = cos(*(fS + 1) * PI / 180);
   XH = *meas * cosb + *(meas+1) * sinb * sina + *(meas+2) * sinb * cosa;
   YH = *(meas+1) * cosa + *(meas+2) * sina;
   return atan2(-YH, XH) * 180 / PI;
}



int main(void)
{
    float meas[9] = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1};
    float sensorErrors[2] = {1., 1.};
    float yaw;
    ParticleFilter<5> filter;
    filter.predict(meas + 3);
    filter.update(getPitchRoll(meas), sensorErrors);
    filter.resample();
    filter.estimate();
    yaw = getYaw(filter.filterState, meas + 6);
    std::cout << "\n" << filter.filterState[0] << "\n" << filter.filterState[1] << "\n";

    return 0;
}