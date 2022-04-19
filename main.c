#include <stdio.h>
#include "inertial.h"

int main() {
    INS_data_t *data;
    INS_state_t *ref;
    INS_result_t state = {91620.0, 23.1373950708 * D2R, 113.3713651222 * D2R, 2.175, 0.0, 0.0, 0.0,
                          0.0107951084511778 * D2R,
                          -2.14251290749072 * D2R, -75.7498049314083 * D2R};
    INS_result_t last = state;
    int i;
    FILE *fout;

//    fout= fopen("ref.txt","w");
//    size_t nref = INS_decode_binary("../../ins/Data1_PureIns.bin", (void *) &ref);
//    for(i=0;i<nref/sizeof (INS_state_t);++i){
//        INS_state(fout,ref+i);
//    }
//    fclose(fout);

    size_t cnt = INS_decode_binary("../../ins/Data1.bin", (void *) &data);
    if (cnt % sizeof(INS_data_t) != 0) {
        printf("data not complete, cut");
    }
    cnt /= sizeof(INS_data_t);
    fout = fopen("result.txt", "w");
    for (i = 0; i < cnt; ++i) {
        if (data[i].sec - state.sec <= 0)
            continue;
        INS_update_param(&data[i], &state,&last);
        INS_update_att(&data[i], &state,&last);
        INS_update_vel(&data[i], &state,&last);
        INS_update_pos(&data[i], &state,&last);
        INS_output(fout, &state);
        last=state;
    }
    fclose(fout);
    free(data);
//    free(ref);
    return 0;
}
