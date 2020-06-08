#include <vector>

namespace etap {
    static const std::vector<double> x_values = {0.01, 0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275, 0.3, 0.325, 0.35, 0.375, 0.4, 0.425, 0.45, 0.475, 0.5, 0.525, 0.55, 0.575, 0.588, 0.6, 0.613, 0.625, 0.63, 0.635, 0.64, 0.645, 0.648, 0.65, 0.653, 0.655, 0.658, 0.66, 0.665, 0.667, 0.668, 0.669, 0.67, 0.671, 0.672, 0.673, 0.674, 0.675, 0.68, 0.685, 0.69, 0.695, 0.7, 0.71, 0.725, 0.75, 0.775, 0.8, 0.825, 0.85, 0.875, 0.9, 0.925, 0.95, 0.975, 0.99};
    static const std::vector<double> y_values = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99};
    static const std::vector<double> corr_values = {-6.6854651, -8.636188, -10.553396000000001, -11.374085, -12.035285, -12.585394, -13.072593, -13.512858000000001, -13.91918, -14.299980000000001, -14.66383, -15.01346, -15.355730000000001, -15.69254, -16.02775, -16.36427, -16.701929999999997, -17.041970000000003, -17.3814, -17.71247, -18.017039999999998, -18.26344, -18.39533, -18.34887, -18.23955, -18.09376, -17.90879, -17.75951, -17.719060000000002, -17.704359999999998, -17.72465, -17.79573, -17.86907, -17.93356, -18.04795, -18.12953, -18.210019999999997, -18.162760000000002, -16.693459999999998, -15.170720000000001, -14.30027, -13.472220000000002, -12.768509999999997, -12.235729999999998, -11.88061, -11.680159999999997, -11.59934, -11.602670000000002, -12.099450000000001, -12.631689999999999, -13.01419, -13.28166, -13.472450000000002, -13.72533, -13.954849999999999, -14.221150000000002, -14.460490000000002, -14.695779999999997, -14.92555, -15.153349999999998, -15.396320000000001, -15.706850000000001, -16.25329, -17.644399999999997, -23.1928, -41.6496, -6.6854452, -8.636177, -10.553467, -11.374283, -12.035622, -12.585863, -13.073175999999998, -13.513528, -13.91991, -14.300740000000001, -14.66458, -15.01418, -15.35638, -15.69311, -16.028200000000002, -16.36458, -16.70208, -17.04193, -17.38116, -17.71202, -18.01636, -18.26251, -18.39415, -18.347430000000003, -18.23797, -18.09207, -17.90699, -17.75763, -17.717160000000003, -17.702469999999998, -17.72279, -17.79395, -17.86736, -17.931900000000002, -18.046390000000002, -18.128049999999998, -18.208679999999998, -18.161510000000003, -16.6924, -15.169700000000002, -14.299259999999999, -13.47122, -12.767509999999998, -12.234719999999998, -11.87958, -11.679109999999998, -11.59826, -11.601560000000001, -12.098120000000002, -12.63012, -13.01238, -13.27965, -13.470260000000001, -13.72287, -13.952059999999998, -14.21795, -14.456930000000002, -14.691859999999998, -14.92124, -15.148539999999999, -15.39075, -15.69987, -16.24299, -17.6241, -23.127299999999998, -41.4156, -6.6853185, -8.6360281, -10.553554, -11.374770999999999, -12.036558, -12.587223, -13.074898, -13.515531000000001, -13.92211, -14.303030000000001, -14.66688, -15.01639, -15.35842, -15.69488, -16.02965, -16.36562, -16.702659999999998, -17.041980000000002, -17.38063, -17.71086, -18.01453, -18.25997, -18.39087, -18.343390000000003, -18.23354, -18.08729, -17.90188, -17.752299999999998, -17.7118, -17.697119999999998, -17.71753, -17.78889, -17.86249, -17.92719, -18.041980000000002, -18.12387, -18.20487, -18.157960000000003, -16.68938, -15.166790000000002, -14.29638, -13.468350000000001, -12.764629999999997, -12.231799999999998, -11.87662, -11.676069999999998, -11.59514, -11.598340000000002, -12.094260000000002, -12.62553, -13.00713, -13.2738, -13.463920000000002, -13.71571, -13.943979999999998, -14.20865, -14.44654, -14.680389999999997, -14.908579999999999, -15.134389999999998, -15.37437, -15.679350000000001, -16.21282, -17.564999999999998, -22.9379, -40.742599999999996, -6.684514, -8.6349428, -10.552956, -11.375052, -12.037769, -12.589282, -13.077665, -13.518845, -13.9258, -14.306930000000001, -14.67083, -15.02024, -15.36203, -15.69811, -16.03238, -16.36774, -16.70405, -17.04256, -17.380300000000002, -17.70954, -18.01215, -18.25647, -18.3862, -18.33752, -18.22706, -18.08026, -17.89432, -17.74439, -17.70382, -17.689159999999998, -17.7097, -17.78136, -17.85525, -17.92019, -18.035410000000002, -18.117639999999998, -18.199189999999998, -18.152650000000005, -16.68481, -15.162370000000001, -14.291989999999998, -13.463960000000002, -12.760199999999998, -12.22731, -11.87203, -11.671369999999998, -11.59029, -11.593330000000002, -12.088220000000002, -12.61835, -12.998899999999999, -13.26466, -13.454000000000002, -13.70452, -13.931309999999998, -14.19401, -14.43012, -14.662169999999998, -14.888349999999999, -15.111679999999998, -15.347990000000001, -15.64633, -16.164649999999998, -17.471999999999998, -22.6437, -39.7006, -6.686081, -8.6366258, -10.5449958, -11.371694999999999, -12.037816, -12.591375000000001, -13.081200999999998, -13.523410000000002, -13.93106, -14.31262, -14.67669, -15.02605, -15.36758, -15.70322, -16.03688, -16.37148, -16.706889999999998, -17.04436, -17.380950000000002, -17.70892, -18.01017, -18.253040000000002, -18.381240000000002, -18.33099, -18.21972, -18.0722, -17.88555, -17.735159999999997, -17.6945, -17.67984, -17.70053, -17.77253, -17.846760000000003, -17.91198, -18.027710000000003, -18.110319999999998, -18.19249, -18.146370000000005, -16.679299999999998, -15.156970000000001, -14.286589999999999, -13.458530000000001, -12.754689999999998, -12.221689999999999, -11.86627, -11.665419999999997, -11.58414, -11.586940000000002, -12.080440000000001, -12.60911, -12.98831, -13.25291, -13.441250000000002, -13.69014, -13.915009999999999, -14.175040000000001, -14.408660000000001, -14.638149999999998, -14.86148, -15.081269999999998, -15.31253, -15.602010000000002, -16.10078, -17.3513, -22.2708, -38.3906, -6.6876429, -8.641161, -10.561705, -11.383633, -12.041812, -12.586044000000001, -13.083925999999998, -13.528870000000001, -13.938, -14.320400000000001, -14.6849, -15.03435, -15.37569, -15.71091, -16.04392, -16.37769, -16.712079999999997, -17.048370000000002, -17.38363, -17.71013, -18.00978, -18.25092, -18.37731, -18.32518, -18.21293, -18.064539999999997, -17.87705, -17.72607, -17.685280000000002, -17.6706, -17.69143, -17.76378, -17.838330000000003, -17.90383, -18.02006, -18.10304, -18.185789999999997, -18.140030000000003, -16.67353, -15.151190000000001, -14.28075, -13.452580000000001, -12.748599999999998, -12.215409999999999, -11.859760000000001, -11.658669999999997, -11.5771, -11.579600000000001, -12.071380000000001, -12.5983, -12.97594, -13.23921, -13.426400000000003, -13.67339, -13.89595, -14.152650000000001, -14.383020000000002, -14.609079999999999, -14.82857, -15.043666999999997, -15.268392, -15.546985000000001, -16.02283, -17.209, -21.846, -36.9156, -6.688409399999999, -8.643443, -10.566817, -11.392661, -12.057851, -12.61032, -13.096768999999998, -13.525150000000002, -13.94601, -14.330760000000001, -14.69627, -15.046109999999999, -15.38744, -15.72233, -16.05474, -16.38766, -16.720989999999997, -17.05602, -17.38983, -17.7147, -18.012539999999998, -18.251730000000002, -18.37603, -18.321720000000003, -18.20834, -18.05894, -17.870449999999998, -17.71877, -17.67781, -17.663059999999998, -17.683970000000002, -17.75659, -17.83142, -17.89714, -18.01376, -18.09702, -18.180179999999996, -18.13462, -16.66819, -15.145620000000001, -14.274999999999999, -13.446600000000002, -12.742359999999998, -12.208879999999999, -11.8529, -11.651449999999999, -11.5695, -11.571606000000001, -12.061299000000002, -12.586229999999999, -12.962159999999999, -13.22397, -13.409910000000002, -13.65479, -13.874669999999998, -14.127270000000001, -14.353470000000002, -14.575019999999999, -14.789385, -14.998282999999997, -15.214712200000001, -15.480289, -15.930499999999999, -17.048199999999998, -21.3915, -35.3686, -6.6891174, -8.645245, -10.570587, -11.398678, -12.066575, -12.622622, -13.114819999999998, -13.558650000000002, -13.96571, -14.340470000000002, -14.71216, -15.06321, -15.40487, -15.73963, -16.07156, -16.40369, -16.73597, -17.06969, -17.40193, -17.725, -18.02083, -18.25779, -18.37967, -18.32282, -18.208099999999998, -18.057499999999997, -17.86782, -17.715269999999997, -17.67406, -17.65917, -17.68006, -17.75282, -17.82779, -17.89363, -18.01042, -18.09378, -18.176969999999997, -18.131330000000002, -16.66405, -15.140850000000002, -14.269849999999998, -13.441040000000001, -12.736349999999998, -12.202389999999998, -11.84591, -11.643940999999998, -11.561464, -11.563035000000001, -12.050126, -12.572794, -12.94685, -13.2071, -13.391690000000002, -13.63425, -13.850999999999999, -14.098450000000001, -14.319110000000002, -14.534509999999997, -14.741832, -14.942321099999997, -15.147926000000002, -15.397684000000002, -15.819439999999998, -16.86687, -20.9199, -33.8166, -6.690008399999999, -8.647361, -10.575039, -11.405438, -12.075775, -12.634500000000001, -13.12979, -13.577470000000002, -13.99022, -14.37583, -14.74119, -15.08954, -15.43267, -15.76767, -16.09934, -16.43082, -16.762069999999998, -17.0944, -17.42492, -17.74592, -18.03934, -18.27358, -18.39246, -18.3324, -18.21598, -18.06385, -17.87264, -17.718919999999997, -17.67734, -17.66216, -17.682850000000002, -17.75547, -17.830350000000003, -17.8961, -18.01266, -18.09575, -18.17831, -18.132030000000004, -16.662219999999998, -15.137650000000002, -14.265899999999998, -13.43631, -12.730829999999997, -12.196049999999998, -11.838757000000001, -11.635973999999997, -11.552695, -11.553479000000001, -12.037039000000002, -12.556937999999999, -12.928848, -13.18736, -13.370450000000002, -13.61032, -13.823149999999998, -14.063550000000001, -14.276230000000002, -14.482533999999998, -14.679378, -14.867488999999997, -15.057780000000001, -15.286830000000002, -15.675540999999999, -16.65112, -20.428, -32.2936, -6.6917696, -8.651693, -10.583596, -11.418156, -12.09261, -12.65546, -13.154959999999999, -13.607020000000002, -14.02448, -14.41544, -14.78794, -15.14376, -15.48659, -15.82213, -16.15403, -16.485129999999998, -16.815389999999997, -17.1461, -17.47437, -17.79251, -18.08246, -18.31265, -18.42695, -18.36191, -18.24284, -18.08831, -17.894669999999998, -17.73902, -17.69675, -17.68094, -17.70102, -17.7729, -17.847140000000003, -17.91232, -18.0277, -18.10971, -18.19005, -18.141820000000003, -16.66546, -15.137720000000002, -14.26433, -13.433090000000002, -12.725949999999997, -12.18954, -11.830641, -11.626307999999998, -11.54154, -11.540905000000002, -12.018541, -12.534308999999999, -12.903297, -13.159546, -13.340676000000002, -13.57681, -13.783649999999998, -14.012236000000001, -14.210848000000002, -14.400777999999999, -14.57866, -14.744609999999998, -14.90844, -15.10453, -15.44803, -16.344829999999998, -19.8655, -30.7646, -6.700015, -8.671114, -10.62134, -11.47378, -12.16556, -12.74526, -13.26113, -13.729080000000002, -14.16191, -14.567720000000001, -14.95456, -15.324349999999999, -15.68272, -16.02938, -16.36417, -16.69636, -17.025329999999997, -17.35238, -17.67465, -17.98451, -18.26396, -18.48153, -18.58133, -18.50036, -18.37277, -18.2105, -18.00899, -17.846899999999998, -17.80216, -17.78387, -17.8011, -17.86901, -17.939760000000003, -18.00186, -18.11104, -18.18754, -18.25705, -18.199700000000004, -16.69426, -15.153240000000002, -14.273119999999999, -13.43518, -12.721459999999997, -12.178651999999998, -11.81361, -11.603422999999998, -11.513132, -11.507331, -11.964613540000002, -12.4676694, -12.8285711, -13.0789332, -13.254953500000003, -13.480701, -13.668861999999999, -13.857431, -14.006400000000001, -14.137609999999999, -14.24709, -14.333609999999998, -14.40501, -14.493210000000001, -14.7105, -15.454709999999999, -18.6752, -28.4346};
}