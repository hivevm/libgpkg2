#ifndef GPKG_I18N_H
#define GPKG_I18N_H

#include <stddef.h>

typedef struct i18n_locale i18n_locale_t;

double i18n_strtod(const char *nptr, char **endptr, i18n_locale_t *locale);

/**
 * snprintf in the given locale, independent of the process's numeric locale. Unlike
 * sqlite3_snprintf, whose floating point conversions cap the significant digits, this formats
 * with the full precision the C library provides.
 */
int i18n_snprintf(char *out, size_t length, i18n_locale_t *locale, const char *format, ...);

i18n_locale_t *i18n_locale_init(const char *locale_name);

void i18n_locale_destroy(i18n_locale_t *locale);

#endif
