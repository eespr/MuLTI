/*

  matlab mex interface to geopsy dispersion calculation

  * usage:

    out = gpdc(T, Vp, Vs, d)
    out = gpdc(T, Vp, Vs, d, 'nSamples', n, 'minRange', mn, 'maxRange', mx)
    out = gpdc(T, Vp, Vs, d, 'fV', f)

  where required input:

    T  = Thickness(m)    : double vector
    Vp = Vp (m/s)        : double vector
    Vs = Vs (m/s)        : double vector
    d  = density (kg/m3) : double vector

  optional input for frequency values. one or more of these options can be set:

    n  = number of samples : scalar integer (but in matlab the type would actually be 'double'...)
    mn = minimum range     : scalar double
    mx = maximum range     : scalar double

  or, the frequency values can be set form a vector of doubles:

    f = frequency vector : double vector

  if 'fV' is provided 'nSamples', 'minRange' and 'maxRange' values will be ignored.

  example using fV:

    out = gpdc(T, Vp, Vs, d, 'fV', [10, 20, 30]); 

  which should be equivalent to:

    gpdc -R 5 -s frequency -n 3 -min 10.0 -max 30.0 test.model

  * output:

    out = 2d matrix

      col1 : x values
      col2 : y values for mode 0
      col3 : y values for mode 1
      ...

  * example:

    T  = [7.5, 25, 0];
    Vp = [500, 1350, 2000];
    Vs = [200, 210, 1000];
    d  = [1700, 1900, 2500];
    out = gpdc(T, Vp, Vs, d);
    plot(out(:, 1), out(:, 2:end))

  which should be equivalent to (these command line options match our defaults):

    gpdc -R 5 -s frequency -n 150 -min 1.0 -max 150.0 test.model | figue -c

*/

/* define name of mex function: */
#define MEXNAME "gpdc"

/* include sstream and string for error messaging: */
#include <sstream>
/* include matlab mex + matrix headers: */
#include "mex.h"
#include "matrix.h"
/* geopsy headers: */
#include <QGpCoreWave.h>


/* error messaging function: */
const char *errorMessage(int msgType,
                         const char *errorMessage) {

  /* stream for error message: */
  std::stringstream errStr;
  /* if type is 0, this is the error id: */
  if (msgType == 0) {
    errStr << "MATLAB:" << MEXNAME << ":" << errorMessage;
  } else {
    /* this is the error message: */
    errStr << MEXNAME << " : " << errorMessage;
  }
  /* to string ... : */
  std::string errMsg = errStr.str();
  /* return message as char: */
  return errMsg.c_str();

}


/* function to process model: */
void processModel(int nLayers,
                  const mxArray *T,
                  const mxArray *Vp,
                  const mxArray *Vs,
                  const mxArray *d,
                  const mxArray *nSamplesIn,
                  const mxArray *minRangeIn,
                  const mxArray *maxRangeIn,
                  const mxArray *fV,
                  mxArray *plhs[]) {

  /* default values: */
  int nRayleigh      = 5;      /* default : 1    */
  int nLove          = 0;
  bool groupSlowness = false;
  int nSamples       = 150;    /* default : 100  */
  double minRange    = 1.0;    /* default : 0.2  */
  double maxRange    = 150.0;  /* default : 20.0 */
  int vNSamples      = 100;
  double vMinRange   = 100.0;
  double vMaxRange   = 3000.0;
  bool oneMode       = false;
  bool force         = false;

  /*
    mode:
      O : GridMode
      1 : CurveMode
  */
  enum AppMode {GridMode = 0, CurveMode = 1};
  AppMode mode = CurveMode;

  /*
    samplingType:
      0 : LinearScale
      1 : LogScale
      2 : InversedScale
      4 : Interpole
      8 : Function
  */
  SamplingOption samplingType = LinearScale; /* default : LogScale */

  /* input: */
  double *T_data  = (double *) mxGetData(T);
  double *Vp_data = (double *) mxGetData(Vp);
  double *Vs_data = (double *) mxGetData(Vs);
  double *d_data  = (double *) mxGetData(d);
  std::vector<double> fV_vector;

  /* variables: */
  int i, j;
  QVector<double> x;

  /* output: */
  Curve<Point2D> curveOut;
  int curveOutCount;
  double *mOut = NULL;

  /* check inputs. if nSamples is provided ... : */
  if (nSamplesIn) {
    /* get supplied value: */
    double *nSamples_data = (double *) mxGetData(nSamplesIn);
    /* should be an integer: */
    if (*nSamples_data != floor(*nSamples_data)) {
      /* or exit: */
      mexErrMsgIdAndTxt(errorMessage(0, "typeargin"),
        errorMessage(1, "nSamples should be an integer"));
    } else {
      /* use supplied value: */
      nSamples = floor(*nSamples_data);
    }
  }

  /* if minRange is not empty ... : */
  if (minRangeIn) {
    /* get supplied value: */
    double *minRangePtr = (double *) mxGetData(minRangeIn);
    minRange = *minRangePtr;
  }

  /* if maxRange is not empty ... : */
  if (maxRangeIn) {
    /* get supplied value: */
    double *maxRangePtr = (double *) mxGetData(maxRangeIn);
    maxRange = *maxRangePtr;
  }

  /* if fV is not empty ... : */
  if (fV) {
    double *fV_data = (double *) mxGetData(fV);
    /* nSamples is length of fV: */
    nSamples = mxGetN(fV);
    /* set fV_vector to appropriate size: */
    fV_vector.resize(nSamples);
    /* spin through fV_data: */
    for (i = 0; i < nSamples; i++) {
      fV_vector[i] = *fV_data;
      *fV_data++;
    }
    /* sort fV_vector: */
    std::sort(fV_vector.begin(), fV_vector.end());
    /* minRange and maxRange are first and last values: */
    minRange = fV_vector[0];
    maxRange = fV_vector[nSamples - 1];
  }

  /* if no plugincoreapplication instance ... : */
  if (! PluginCoreApplication::instance()) {
    /* set up plugincoreapplication: */
    new PluginCoreApplication;
    /* disable application output (stout / stderr): */
    PluginCoreApplication::instance()->freezeStream(true);
  }

  /* initialise model: */
  LayeredModel model(nLayers);

  /* add data to model: */
  for (i = 0; i < nLayers; i++) {
    /* don't set thickness for final layer: */
    if(i < nLayers-1) {
      /* add T: */
      model.setH(i, *T_data);
      *T_data++;
    }
    /* add Vp, Vs, d: */
    model.setSlowP(i, 1.0 / *Vp_data);
    *Vp_data++;
    model.setSlowS(i, 1.0 / *Vs_data);
    *Vs_data++;
    model.setRho(i, *d_data);
    *d_data++;
    /* non mandatory quality factors, presume 0: */
    model.setQp(i, 0.0);
    model.setQs(i, 0.0);
  }

  /* initialise model calculation ... : */
  model.initCalculation();

  /* compute common sampling scale: */
  Curve<Point1D> curve;
  curve.line(minRange, 0.0, maxRange, 0.0);
  curve.resample(nSamples, minRange, maxRange, samplingType | Function);

  /* if fV has been supplied ... : */
  if (fV) {
    /* spin through the frequency vector: */
    for (i = 0; i < nSamples; i++) {
      /* set the curve x value: */
      curve[i].setX(fV_vector[i]);
    }
  }

  /* convert to angular frequency: */
  curve.xMultiply(2*M_PI);
  x = curve.xVector();

  /* mode switching: */
  switch(mode) {

    /* curve mode (1): */
    case CurveMode:

      /* check number of rayleigh modes: */
      if (nRayleigh > 0) {
        /* rayleigh the model ... : */
        Rayleigh rayleigh(&model);
        /* initialise dispersion: */
        Dispersion dispersion(nRayleigh, &x);
        /* if dispersion can be calculated ... : */
        if (dispersion.calculate(&rayleigh, 0)) {

          /* if using groupSlowness: */
          if (groupSlowness) {
            dispersion.setGroupSlowness();
          }

          /* send output curves back to matlab ... init output matrix : x + nRaylegh*y : */
          plhs[0] = mxCreateDoubleMatrix(nSamples, (nRayleigh + 1), mxREAL);
          /* get pointer for output: */
          mOut = mxGetPr(plhs[0]);
          /* for each rayleigh mode: */
          for (i = 0; i < nRayleigh; i++) {
            /* get dispersion curve values ... : */
            curveOut = dispersion.curve(i);
            /* get curve value count: */
            curveOutCount = curveOut.count();
            /* if this is mode 0: */
            if (i == 0) {
              /* set x values, 1 -> nSamples: */
              for (j = 0; j < curveOutCount; j++) {
                const Point2D p = curveOut.at(j);
                *mOut = p.x();
                *mOut++;
              }
            }
            /* set any null values ... : */
            for (j = 0; j < (nSamples - curveOutCount); j++) {
              *mOut = mxGetNaN();
              *mOut++;
            }
            /* loop through curve, and set y values: */
            for (j = 0; j < curveOutCount; j++) {
              const Point2D p = curveOut.at(j);
              *mOut = p.y();
              *mOut++;
            }
          }
          /* end matlab output. */

        } else if(force) {
          /* cannot compute dispersion curves ... : */
          mexErrMsgIdAndTxt(errorMessage(0, "computedispersioncurves"),
                            errorMessage(1, "cannot compute dispersion curves"));
        } else {
          /* dispersion calculate failed ... : */
          mexErrMsgIdAndTxt(errorMessage(0, "dispersioncalculate"),
                            errorMessage(1, "dispersion calculate failed"));
        }

      } else {
        /* rayleigh number not greater than 0, exit: */
        mexErrMsgIdAndTxt(errorMessage(0, "typeargin"),
                          errorMessage(0, "rayleigh should be > 0 for curve mode"));
      }

      /* curve mode done: */
      break;

    /* default is exit: */
    default:

      /* valid mode not specified: */
      mexErrMsgIdAndTxt(errorMessage(0, "typeargin"),
                        errorMessage(1, "mode specified is not valid"));
      break;

  }

  /* return: */
  return;

}


/* main gateway mexFunction: */
void mexFunction(int nlhs, mxArray* plhs[],
                 int nrhs, const mxArray* prhs[]) {

  /* input variables: */
  const mxArray *T        = NULL;
  const mxArray *Vp       = NULL;
  const mxArray *Vs       = NULL;
  const mxArray *d        = NULL;
  const mxArray *nSamples = NULL;
  const mxArray *minRange = NULL;
  const mxArray *maxRange = NULL;
  const mxArray *fV       = NULL;
  int nLayers;

  /* variables: */
  int i;
  std::stringstream errStr, s;

  /* check number of output arguments is either 0 or 1: */
  if ((nlhs != 0) &&
      (nlhs != 1)) {
    /* or exit: */
    mexErrMsgIdAndTxt(errorMessage(0, "nargout"),
                      errorMessage(1, "one output argument expected"));
  }

  /* check number input of arguments is 4, 6, 8, 10 or 12: */
  if ((nrhs != 4) &&
      (nrhs != 6) &&
      (nrhs != 8) &&
      (nrhs != 10) &&
      (nrhs != 12)) {
    /* or exit: */
    mexErrMsgIdAndTxt(errorMessage(0, "nargin"),
                      errorMessage(1, "incorrect number of input arguments"));
  }

  /* check input arguments ... T: */
  if (!mxIsDouble(prhs[0]) ||
     mxIsScalar(prhs[0]) ||
     mxIsComplex(prhs[0])) {
    /* or exit: */
    mexErrMsgIdAndTxt(errorMessage(0, "typeargin"),
                      errorMessage(1, "first input should be double vector : Thickness (m)"));
  } else {
    /* use supplied value: */
    T = prhs[0];
  }

  /* check input arguments ... Vp: */
  if (!mxIsDouble(prhs[1]) ||
     mxIsScalar(prhs[1]) ||
     mxIsComplex(prhs[1])) {
    /* or exit: */
    mexErrMsgIdAndTxt(errorMessage(0, "typeargin"),
                      errorMessage(1, "second input should be double vector : Vp (m/s)"));
  } else {
    /* use supplied value: */
    Vp = prhs[1];
  }

  /* check input arguments ... Vs: */
  if (!mxIsDouble(prhs[2]) ||
     mxIsScalar(prhs[2]) ||
     mxIsComplex(prhs[2])) {
    /* or exit: */
    mexErrMsgIdAndTxt(errorMessage(0, "typeargin"),
                      errorMessage(1, "third input should be double vector : Vs (m/s)"));
  } else {
    /* use supplied value: */
    Vs = prhs[2];
  }

  /* check input arguments ... d: */
  if (!mxIsDouble(prhs[3]) ||
     mxIsScalar(prhs[3]) ||
     mxIsComplex(prhs[3])) {
    /* or exit: */
    mexErrMsgIdAndTxt(errorMessage(0, "typeargin"),
                      errorMessage(1, "fourth input should be double vector : density (kg/m3)"));
  } else {
    /* use supplied value: */
    d = prhs[3];
  }

  /*
    arguments > 3 should be optional key / value pairs, i.e,:

      'nSamples', 100
      'minRange', 1.0
      'maxRange', 150.0
  */

  /* while 3 < i < number of arguments ... : */
  for (i = 4; i < nrhs; i += 2) {

    /* check this argument is a char: */
    if (mxIsChar(prhs[i])) {
      /* check length of argument: */
      mwSize argLen = mxGetN(prhs[i]) + 1;
      /* if less than 16 ... : */
      if (argLen < 16) {
        /* get argName: */
        char *argName = new char[16];
        mxGetString(prhs[i], argName, argLen);

        /* if nSamples ... : */
        if (strcmp(argName, "nSamples") == 0) {
          /* argValue is i+1, and should be an integer: */
          if (!mxIsDouble(prhs[i + 1]) ||
              !mxIsScalar(prhs[i + 1]) ||
              mxIsComplex(prhs[i + 1])) {
            /* or exit: */
            mexErrMsgIdAndTxt(errorMessage(0, "typeargin"),
                              errorMessage(1, "input nSamples should be an integer"));
          } else {
            /* use supplied value: */
            nSamples = prhs[i + 1];
          }

        /* if minRange: */
        } else if (strcmp(argName, "minRange") == 0) {
          /* argValue is i+1, and should be a double: */
          if (!mxIsDouble(prhs[i + 1]) ||
              !mxIsScalar(prhs[i + 1]) ||
              mxIsComplex(prhs[i + 1])) {
            /* or exit: */
            mexErrMsgIdAndTxt(errorMessage(0, "typeargin"),
                              errorMessage(1, "input minRange should be a double"));
          } else {
            /* use supplied value: */
            minRange = prhs[i + 1];
          }

        /* if maxRange: */
        } else if (strcmp(argName, "maxRange") == 0) {
          /* argValue is i+1, and should be a double: */
          if (!mxIsDouble(prhs[i + 1]) ||
              !mxIsScalar(prhs[i + 1]) ||
              mxIsComplex(prhs[i + 1])) {
            /* or exit: */
            mexErrMsgIdAndTxt(errorMessage(0, "typeargin"),
                              errorMessage(1, "input maxRange should be a double"));
          } else {
            /* use supplied value: */
            maxRange = prhs[i + 1];
          }

        /* if fV (frequency vector): */
        } else if (strcmp(argName, "fV") == 0) {
          /* argValue is i+1, and should be a double vector: */
          if (!mxIsDouble(prhs[i + 1]) ||
              mxIsScalar(prhs[i + 1]) ||
              mxIsComplex(prhs[i + 1])) {
            /* or exit: */
            mexErrMsgIdAndTxt(errorMessage(0, "typeargin"),
                              errorMessage(1, "input fV should be double vector"));
          } else {
            /* use supplied value: */
            fV = prhs[i + 1];
            /* fV nullifies nSamples, minRange and maxRange: */
            nSamples = NULL;
            minRange = NULL;
            maxRange = NULL;
          }

        } else {
          /* invalid argument name ... convert i to string: */
          s << i;
          /* create message: */
          errStr << "invalid argument name '" << argName << "' at position " << s.str();
          /* error and exit ... : */
          mexErrMsgIdAndTxt(errorMessage(0, "nameargin"),
                            errorMessage(1, errStr.str().c_str()));
        }

      } else {
        /* invalid argument length (char > 16) ... convert i to string: */
          s << i;
        /* create message: */
        errStr << "invalid argument length at position " << s.str();
        /* error and exit ... : */
        mexErrMsgIdAndTxt(errorMessage(0, "typeargin"),
                          errorMessage(1, errStr.str().c_str()));
      }

    } else {
      /* invalid argument type ... convert i to string: */
      s << i;
      /* create message: */
      errStr << "invalid input argument type at position " << s.str();
      /* error and exit ... : */
      mexErrMsgIdAndTxt(errorMessage(0, "typeargin"),
                        errorMessage(1, errStr.str().c_str()));
    }

  }
  /* end argument > 3 */

  /* check inputs are of the same length: */
  if (mxGetN(T)  != mxGetN(Vp) ||
      mxGetN(Vp) != mxGetN(Vs) ||
      mxGetN(Vs) != mxGetN(d)) {
    /* or exit: */
    mexErrMsgIdAndTxt(errorMessage(0, "lenargin"),
                      errorMessage(1, "inputs T, Vp, Vs and d should be of equal length"));
  }

  /* number of layers: */
  nLayers = mxGetN(T);
  /* process model - function sends output back to matlab: */
  processModel(nLayers, T, Vp, Vs, d, nSamples, minRange, maxRange, fV, plhs);
  /* return: */
  return;

}

