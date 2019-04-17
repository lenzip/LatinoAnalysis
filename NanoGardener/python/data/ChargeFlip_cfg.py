ChargeFlip={

# --------------------------- Full2017 ---------------------------------

  'Full2017v2'  :  {
                     #                  |eta| bin    pT bin         DATA       ERR        SysErrDown SysErrUp    MC          ERR       SF      ERR      SysErrDown SysErrUp
                     'FlipProba' :   [
                                      [ 0.0 , 0.5 ,  13.0 , 20.0  , 1.02E-3  , 0.14E-3  , 0.78E-3   , 0.66E-3  , 2.69E-5   , 0.69E-5  , 38.0  , 11.0  , 29        , 24       ],
                                      [ 0.5 , 1.0 ,  13.0 , 20.0  , 1.01E-3  , 0.13E-3  , 0.53E-3   , 0.65E-3  , 0.912E-4  , 0.090E-4 , 11.0  , 1.8   , 7.1       , 5.9      ],
                                      [ 1.0 , 1.5 ,  13.0 , 20.0  , 1.40E-3  , 0.15E-3  , 0.66E-3   , 0.65E-3  , 3.09E-4   , 0.18E-4  , 4.54  , 0.55  , 2.2       , 2.1      ],
                                      [ 1.5 , 2.0 ,  13.0 , 20.0  , 2.66E-3  , 0.18E-3  , 1.1E-3    , 1.1E-3   , 9.88E-4   , 0.31E-4  , 2.70  , 0.20  , 1.1       , 1.1      ],
                                      [ 2.0 , 2.5 ,  13.0 , 20.0  , 2.90E-3  , 0.22E-3  , 1.2E-3    , 1.2E-3   , 1.488E-3  , 0.046E-3 , 1.64  , 0.16  , 0.45      , 1.1      ],
                                      [ 0.0 , 0.5 ,  20.0 , 200.0 , 4.85E-5  , 0.49E-5  , 1.1E-5    , 1.2E-5   , 2.69E-5   , 0.69E-5  , 1.80  , 0.50  , 0.38      , 0.42     ],
                                      [ 0.5 , 1.0 ,  20.0 , 200.0 , 1.306E-4 , 0.064E-4 , 0.28E-4   , 0.29E-4  , 0.912E-4  , 0.090E-4 , 1.43  , 0.16  , 0.30      , 0.31     ],
                                      [ 1.0 , 1.5 ,  20.0 , 200.0 , 4.80E-4  , 0.12E-4  , 1.1E-4    , 1.1E-4   , 3.09E-4   , 0.18E-4  , 1.551 , 0.095 , 0.33      , 0.33     ],
                                      [ 1.5 , 2.0 ,  20.0 , 200.0 , 1.522E-3 , 0.021E-3 , 0.33E-3   , 0.32E-3  , 9.88E-4   , 0.31E-4  , 1.540 , 0.053 , 0.33      , 0.32     ],
                                      [ 2.0 , 2.5 ,  20.0 , 200.0 , 2.295E-3 , 0.033E-3 , 0.49E-3   , 0.47E-3  , 1.488E-3  , 0.046E-3 , 1.542 , 0.052 , 0.32      , 0.32     ],
                                     ],

#                    #                  |eta| bin    pT bin         DATA       ERR        MC          ERR        SF      ERR
#                    'FlipProba' : [
#                                     [ 0.0 , 0.5 ,  13.0 , 20.0  , 1.02E-3  , 0.14E-3  , 2.69E-5   , 0.69E-5  , 38.0  , 11.0  ],
#                                     [ 0.5 , 1.0 ,  13.0 , 20.0  , 1.01E-3  , 0.13E-3  , 0.912E-4  , 0.090E-4 , 11.0  , 1.8   ],
#                                     [ 1.0 , 1.5 ,  13.0 , 20.0  , 1.40E-3  , 0.15E-3  , 3.09E-4   , 0.18E-4  , 4.54  , 0.55  ],
#                                     [ 1.5 , 2.0 ,  13.0 , 20.0  , 2.66E-3  , 0.18E-3  , 9.88E-4   , 0.31E-4  , 2.70  , 0.20  ],
#                                     [ 2.0 , 2.5 ,  13.0 , 20.0  , 2.90E-3  , 0.22E-3  , 1.488E-3  , 0.046E-3 , 1.64  , 0.16  ],
#                                     [ 0.0 , 0.5 ,  20.0 , 200.0 , 4.85E-5  , 0.49E-5  , 2.69E-5   , 0.69E-5  , 1.80  , 0.50  ],
#                                     [ 0.5 , 1.0 ,  20.0 , 200.0 , 1.306E-4 , 0.064E-4 , 0.912E-4  , 0.090E-4 , 1.43  , 0.16  ],
#                                     [ 1.0 , 1.5 ,  20.0 , 200.0 , 4.80E-4  , 0.12E-4  , 3.09E-4   , 0.18E-4  , 1.551 , 0.095 ],
#                                     [ 1.5 , 2.0 ,  20.0 , 200.0 , 1.522E-3 , 0.021E-3 , 9.88E-4   , 0.31E-4  , 1.540 , 0.053 ],
#                                     [ 2.0 , 2.5 ,  20.0 , 200.0 , 2.295E-3 , 0.033E-3 , 1.488E-3  , 0.046E-3 , 1.542 , 0.052 ]
#                                  ],
                      #
                      'SSOSMC' : { 'DY'  : {
                                              '0j'  : [ 0.001052 , 0.000066 ],
                                              '1j'  : [ 0.00120  , 0.00014  ],
                                              '2j'  : [ 0.00185  , 0.00034  ],
                                              'vbs' : [ 0.0041   , 0.0024   ],
                                           }, 
                                   'WW'  : {
                                              '0j'  : [ 0.00126  , 0.00025  ],
                                              '1j'  : [ 0.0098   , 0.00032  ],
                                              '2j'  : [ 0.00273  , 0.00092  ],
                                              'vbs' : [ 0.0041   , 0.0024   ],
                                           },
                                   'Top' : {
                                              '0j'  : [ 0.01084 , 0.0013   ],
                                              '1j'  : [ 0.00584 , 0.00035  ],
                                              '2j'  : [ 0.00354 , 0.00020  ],
                                              'vbs' : [ 0.0041  , 0.0024   ],
                                           },
                                 },
                   },

}
