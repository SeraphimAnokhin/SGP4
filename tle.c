#pragma once

#include "stdio.h"
#include "common.h"


int get_elems(tle_set_t *tle) {
    char line1[71], line2[71];
    FILE *f = fopen("tle.txt", "r");
    fgets(line1, 71, f);
    fgets(line2, 71, f);
    fclose(f);

    if ((line1[0] != '1') || (line2[0] != '2') || (line1[1] != ' ') || (line2[1] != ' ')) {
        return 1; // invalid data line; valid one must have line number as the zero characher and a blank as the first
    }
    int sum1 = 0, sum2 = 0, i;
    for (i = 0; i < 68; i++) {
        if ((line1[i] >= '1') && (line1[i] <= '9')) {
            sum1 += line1[i] - '0';
        }
        else if (line1[i] == '-') {
            sum1++;
        }

        if ((line2[i] >= '1') && (line2[i] <= '9')) {
            sum2 += line2[i] - '0';
        }
        else if (line2[i] == '-') {
            sum2++;
        }
    }
    if (((sum1 % 10) != (line1[68] - '0')) || ((sum2 % 10) != (line2[68] - '0'))) {
        return 2; // wrong line checksum
    }

    int year = (line1[18] - '0') * 10 + line1[19] - '0';
    double day = (line1[20] - '0') * 100 + (line1[21] - '0') * 10 + (line1[22] - '0') + (line1[24] - '0') * 1e-1
        + (line1[25] - '0') * 1e-2 + (line1[26] - '0') * 1e-3 + (line1[27] - '0') * 1e-4 + (line1[28] - '0') * 1e-5
        + (line1[29] - '0') * 1e-6 + (line1[30] - '0') * 1e-7 + (line1[31] - '0') * 1e-8;
    tle->epoch = (double)year * 1e3 + day;

    int B_star_exp = ((line1[59] == '-') ? -1 : 1) * (line1[60] - '0');
    tle->B_star = ((line1[53] == '-') ? -1.0 : 1.0) * ((line1[54] - '0') * 1e-1 + (line1[55] - '0') * 1e-2
        + (line1[56] - '0') * 1e-3 + (line1[57] - '0') * 1e-4 + (line1[58] - '0') * 1e-5) * pow(10, B_star_exp);

    tle->i = ((line2[8] == ' ') ? 0.0 : (line2[8] - '0')) * 100 + ((line2[9] == ' ') ? 0.0 : (line2[9] - '0')) * 10
        + line2[10] - '0' + (line2[12] - '0') * 1e-1 + (line2[13] - '0') * 1e-2 + (line2[14] - '0') * 1e-3
        + (line2[15] - '0') * 1e-4;

    tle->node = ((line2[17] == ' ') ? 0.0 : (line2[17] - '0')) * 100 + ((line2[18] == ' ') ? 0.0 : (line2[18] - '0')) * 10
        + line2[19] - '0' + (line2[21] - '0') * 1e-1 + (line2[22] - '0') * 1e-2 + (line2[23] - '0') * 1e-3
        + (line2[24] - '0') * 1e-4;

    tle->e = (line2[26] - '0') * 1e-1 + (line2[27] - '0') * 1e-2 + (line2[28] - '0') * 1e-3 + (line2[29] - '0') * 1e-4
        + (line2[30] - '0') * 1e-5 + (line2[31] - '0') * 1e-6 + (line2[32] - '0') * 1e-7;
    
    tle->omega = ((line2[34] == ' ') ? 0.0 : (line2[34] - '0')) * 100 + ((line2[35] == ' ') ? 0.0 : (line2[35] - '0')) * 10
        + line2[36] - '0' + (line2[38] - '0') * 1e-1 + (line2[39] - '0') * 1e-2 + (line2[40] - '0') * 1e-3
        + (line2[41] - '0') * 1e-4;
    
    tle->M = ((line2[43] == ' ') ? 0.0 : (line2[43] - '0')) * 100 + ((line2[44] == ' ') ? 0.0 : (line2[44] - '0')) * 10
        + line2[45] - '0' + (line2[47] - '0') * 1e-1 + (line2[48] - '0') * 1e-2 + (line2[49] - '0') * 1e-3
        + (line2[50] - '0') * 1e-4;
    
    tle->n = ((line2[52] == ' ') ? 0.0 : (line2[52] - '0')) * 10 + line2[53] - '0' + (line2[55] - '0') * 1e-1
        + (line2[56] - '0') * 1e-2 + (line2[57] - '0') * 1e-3 + (line2[58] - '0') * 1e-4 + (line2[59] - '0') * 1e-5
        + (line2[60] - '0') * 1e-6 + (line2[61] - '0') * 1e-7 + (line2[62] - '0') * 1e-8;
    return 0;
}