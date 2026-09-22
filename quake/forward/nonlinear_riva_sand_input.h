#ifndef NONLINEAR_RIVA_SAND_INPUT_H
#define NONLINEAR_RIVA_SAND_INPUT_H

/* Host-only, row-aware input normalization. Never use Vs_limit as the mesh
 * wave velocity. G0=0 is an input sentinel, not a zero constitutive modulus. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "nonlinear_riva_sand.h"

/* Strict optional scalar reader shared by field and material-point inputs.
 * A missing selector means type 1; an explicit selector is valid for 02 only.
 * Reject duplicate keys and numeric prefixes such as 1.5 or 2junk.
 * Reject retired keys rather than silently selecting the default. */
static inline int riva_read_reversal_type(FILE *fp, const char *key,
    int is_02, int32_t *out)
{
    if (!fp || !key || !out) return -1;
    rewind(fp);
    char line[4096];
    int found=0, selected=1;
    const char *retired_keys[]={"reversalType","rivasand02_reversalType"};
    while (fgets(line,sizeof(line),fp)) {
        if (!strchr(line,'\n') && !feof(fp)) return -1;
        char *comment=strchr(line,'#'); if (comment) *comment='\0';
        char *s=line; while (isspace((unsigned char)*s)) ++s;
        for (const char *retired : retired_keys) {
            const size_t retired_n=strlen(retired);
            if (!strncmp(s,retired,retired_n) &&
                (!s[retired_n] || isspace((unsigned char)s[retired_n]) || s[retired_n]=='='))
                return -1;
        }
        const size_t n=strlen(key);
        if (strncmp(s,key,n) || (s[n] && !isspace((unsigned char)s[n]) && s[n]!='='))
            continue;
        if (found++ || !is_02) return -1;
        s+=n; while (isspace((unsigned char)*s)) ++s;
        if (*s++!='=') return -1;
        while (isspace((unsigned char)*s)) ++s;
        if (*s!='1' && *s!='2') return -1;
        selected=*s++-'0';
        while (isspace((unsigned char)*s)) ++s;
        if (*s) return -1;
    }
    if (ferror(fp)) return -1;
    *out=selected;
    return 0;
}

static inline int riva_read_material_rows(FILE *fp, const char *key,
    int rows, int custom, int with_vs, const riva_parameters_t *p, double *out)
{
    if (!fp || !p || !out || rows<1) return -1;
    rewind(fp);
    char line[4096];
    int found=0, row=0, width=0, automatic=-1;
    const int stride=12+with_vs;
    while (fgets(line,sizeof(line),fp)) {
        if (!strchr(line,'\n') && !feof(fp)) return -1;
        char *comment=strchr(line,'#'); if (comment) *comment='\0';
        char *s=line; while (isspace((unsigned char)*s)) ++s;
        if (!found) {
            const size_t n=strlen(key);
            if (strncmp(s,key,n) || !(isspace((unsigned char)s[n]) || s[n]=='=')) continue;
            s+=n; while (isspace((unsigned char)*s)) ++s;
            if (*s!='=') return -1;
            ++s; found=1;
        }
        while (isspace((unsigned char)*s)) ++s;
        if (!*s) continue;
        if (row==rows) {
            // Another numeric row means the declared row count was wrong.
            char *end=NULL; (void)strtod(s,&end);
            return end==s ? 0 : -1;
        }
        double values[13]; int count=0;
        while (*s) {
            if (count==13) return -1;
            char *end=NULL; const double value=strtod(s,&end);
            if (end==s || !isfinite(value) ||
                (*end && !isspace((unsigned char)*end))) return -1;
            values[count++]=value; s=end;
            while (isspace((unsigned char)*s)) ++s;
        }
        if (custom ? count!=stride :
            (count!=1+with_vs && count!=2+with_vs)) return -1;
        if (width && width!=count) return -1;
        width=count;
        double *v=out+row*stride;
        if (with_vs) v[0]=values[0];
        v+=with_vs;
        const double defaults[12]={0,0,p->M,p->kd,p->h,p->m,p->zeta,
            p->e_max,p->e_min,p->Q,p->R,p->n_G};
        memcpy(v,defaults,sizeof(defaults));
        for (int i=with_vs;i<count;++i) v[i-with_vs]=values[i];
        if (v[0]<0 || v[0]>1 || v[1]<0) return -1;
        const int mode=v[1]==0.0;
        // Never interpolate an automatic sentinel with a positive coefficient.
        if (automatic>=0 && automatic!=mode) return -1;
        automatic=mode;
        ++row;
    }
    return found && row==rows ? 0 : -1;
}
#endif
