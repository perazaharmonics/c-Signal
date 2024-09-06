  //  Currently the DrawTelfmtStatus() and DrawRateGraph() methods are using y0
  // for the instantaneous value and y1 for the running average.
  //  X-Windows might have trouble keeping up if you make nmax bigger than 512,
  // even though some passes will take longer than 512 seconds.
  class Rates
  {
  public:
    inline void init (void) { n=first=0;last=-1; }
    inline void append (const float y0,const float y1,const bool flag)
    {
      if (n<nmax)                       // Are we buffering?
        n++;                            // Yes.
      else                              // Else we're now a FIFO.
        first=++first%nmax;             // So move the first item pointer.
      last=++last%nmax;                 // Always move the last item pointer.
      y[last][0]=y0;                    // Store the y0 value.
      y[last][1]=y1;                    // Store the y1 value.
      flags[last]=flag;                 // Store the flag value.
    }
    static constexpr int nmax{512};     // Maximum number of values to store.
    float y[nmax][2];                   // The sets of y values.
    bool flags[nmax];                   // The flags.
    int n{0};                           // Number of values stored.
    int first{0},last{-1};              // Indices of first and last values.
  } rates;
  class Averages
  {
  public:
    inline void init (void) { n=first=0;last=-1; }
    inline void append (const float rate)
    {
      if (n<nmax)
        n++;
      else
        first=++first%nmax;
      last=++last%nmax;
      y[last]=rate;
    }
    inline float get (void) const
    {
      if (n==0)
        return 0.0f;
      float sum{y[first]};
      for (int i=1,j=first;i<n;i++)
      {
        j=++j%nmax;
        sum+=y[j];
      }
      return sum/float(n);
    }
    static constexpr int nmax{24};
    float y[nmax];                      // The y values.
    int n{0};                           // Number of values stored.
    int first{0},last{-1};              // Indices of first and last values.
  } rateavg;
  // Measure the similarity between two signals using Cross-Correlations
  class CrossCorrelation
  {
  public:
    inline void init (void) {n=first=0;last=-1; }
    // Store pairs of samples for y0 and y1 (unknown) in a circular buffer.
    inline void append (const float y0, const float y1)
    {
      if (n<nmax)                       // Are we buffering?
        n++;                            // Yes. Keep counting buffer items.
      else                              // Else, we're full. We're a FIFO.
        first=++first%max;              // Move the first item pointer.
      last=++last%max;                  // Always move the last item pointer.
      y[last][0]=y0;                    // Store signal 1 sample.
      y[last][1]=y1;                    // Store signal 2 sample.
    }
    // Calculate the Cross-Correlation between the two signals by dividing by 
    //  their standard deviation product and unbiasing by (1/(n-1)).
    // Reference: https://en.wikipedia.org/wiki/Cross-correlation
    // https://www.dsprelated.com/freebooks/mdft/Unbiased_Cross_Correlation.html
    inline float get (void) const
    {
      if (n==0)                           // Do we have any samples?
        return 0.0f;                      // No, so return 0.
      float sum{0.0f};                    // Initialize sum for correlation.
      float meany0{0.0f};                 // Initialize the mean of y0.
      float meany1{0.0f};                 // Initialize the mean of y1.
      float vary0{0.0f};                  // Initialize variance of y0.
      float vary1{0.0f};                  // Initialize variance of y1.
      for (int i=0,j=first;i<n;i++)       // Loop through all available samples..
      {                                   // .. and recollect their values..
        j=(first+i)%nmax;                 // Current circular buffer index.
        meany0+=y[j][0];                  // Accumulate y0's values for its mean.
        meany1+=y[j][1];                  // Accumulate y1's values for its mean.
      }
      // Compute the Expected value.
      meany0/=float(n);                   // Total samples divided by # of samples.
      meany1/=float(n);                   // Total samples divided by # of samples.
      for (int i=0,j=first;i<n;i++)       // Loop through all available values of
      {                                   // the signals and get std. deviations. 
        j=(first+i)%nmax;                 // Current circular buffer index.
        float devy0=(y[j][0]-meany0);     // Deviation of y0 from the mean.
        float devy1=(y[j][1]-meany1);     // Deviation of y1 from the mean.
        sum+=devy0*devy1;                 // Accumulate product of deviations.
        vary0+=devy0*devy0;               // Variance of y0 is its std. dev squared.
        vary1+=devy1*devy1;               // Variance of y1 is its std. dev squared.               
      }
      // Return unbiased and normalized cross-correlation of y0 and y1.
      return sum/(sqrt(vary0*vary1)*(n-1));
    }
    static constexpr int nmax{512};       // Maximum number of values to store.
    float y[nmax][2];                     // Samples of y0, y1.
    int n{0};                             // Number of values in the buffer.
    int first{0},last{-1};                // Indices of values in the FIFO.
    } ccorr;
  class CrossCovariance
  {
  public:
    inline void append (const float y0, const float y1)
    {
      if (n<nmax)                       // Are we buffering?
        n++;                            // Yes. Keep counting buffer items.
      else                              // Else, we're full. We're a FIFO.
        first=++first%max;              // Move the first item pointer.
      last=++last%max;                  // Always move the last item pointer.
      y[last][0]=y0;                    // Store signal 1 sample.
      y[last][1]=y1;                    // Store signal 2 sample.
    }
    // Calculate the unbiased Cross-Covariance between the two signal.
    // Reference: https://en.wikipedia.org/wiki/Cross-correlation
    // https://www.dsprelated.com/freebooks/mdft/Unbiased_Cross_Correlation.html
    inline float get (void) const
    {
      if (n==0)                           // Do we have any samples?
        return 0.0f;                      // No, so return 0.
      float sum{0.0f};                    // Initialize sum for correlation.
      float meany0{0.0f};                 // Initialize the mean of y0.
      float meany1{0.0f};                 // Initialize the mean of y1.
      for (int i=0,j=first;i<n;i++)       // Loop through all available samples..
      {                                   // .. and recollect their values..
        j=(first+i)%nmax;                 // Current circular buffer index.
        meany0+=y[j][0];                  // Accumulate y0's values for its mean.
        meany1+=y[j][1];                  // Accumulate y1's values for its mean.
      }
      // Compute the Expected value.
      meany0/=float(n);                   // Total samples divided by # of samples.
      meany1/=float(n);                   // Total samples divided by # of samples.
      for (int i=0,j=first;i<n;i++)       // Loop through all available values of
      {                                   // the signals and get std. deviations. 
        j=(first+i)%nmax;                 // Current circular buffer index.
        sum+=(y[j][0]-meany0)*(y[j][1]-meany1);// The Cross-Covariance.              
      }
      // Return the unbiased cross-covariance of y0 and y1.
      return sum/float(n-1);
    }
    static constexpr int nmax{512};       // Maximum number of values to store.
    float y[nmax][2];                     // Samples of y0, y1.
    int n{0};                             // Number of values in the buffer.
    int first{0},last{-1};                // Indices of values in the FIFO.
  } ccovar;
