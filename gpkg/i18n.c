#include "i18n.h"

#ifdef GPKG_HAVE_CONFIG_H
#include "config.h"
#endif

#include "sqlite3.h"

#if defined(LOCALE_USE__CREATE_LOCALE)

#include <locale.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>

struct i18n_locale {
  _locale_t locale;
};

double i18n_strtod(const char *nptr, char **endptr, i18n_locale_t *locale) {
  return _strtod_l(nptr, endptr, locale->locale);
}

int i18n_snprintf(char *out, size_t length, i18n_locale_t *locale, const char *format, ...) {
  va_list args;
  va_start(args, format);
  int result = _vsnprintf_l(out, length, format, locale->locale, args);
  va_end(args);
  /* _vsnprintf_l does not terminate the buffer when the output is truncated. */
  if (length > 0) {
    out[length - 1] = 0;
  }
  return result;
}

i18n_locale_t *i18n_locale_init(const char *locale_name) {
  _locale_t locale;
  i18n_locale_t *locale_struct;

  locale_struct = (i18n_locale_t *)sqlite3_malloc(sizeof(i18n_locale_t));
  if (locale_struct == NULL) {
    return NULL;
  }

  locale = _create_locale(LC_ALL, "C");
  if (locale == NULL) {
    sqlite3_free(locale_struct);
    return NULL;
  }

  locale_struct->locale = locale;
  return locale_struct;
}

void i18n_locale_destroy(i18n_locale_t *locale) {
  if (locale != NULL) {
    _free_locale(locale->locale);
    locale->locale = NULL;

    sqlite3_free(locale);
  }
}

#elif defined(LOCALE_USE_NEWLOCALE)

#define _POSIX_C_SOURCE 200809L
#define _GNU_SOURCE

#ifdef HAVE_LOCALE_H
#include <locale.h>
#endif

#ifdef HAVE_XLOCALE_H
#include <xlocale.h>
#endif

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>

struct i18n_locale {
  locale_t locale;
};

double i18n_strtod(const char *nptr, char **endptr, i18n_locale_t *locale) {
  return strtod_l(nptr, endptr, locale->locale);
}

int i18n_snprintf(char *out, size_t length, i18n_locale_t *locale, const char *format, ...) {
  va_list args;
  va_start(args, format);
  /* uselocale only switches the calling thread's locale, so this is thread safe. */
  locale_t previous = uselocale(locale->locale);
  int result = vsnprintf(out, length, format, args);
  uselocale(previous);
  va_end(args);
  return result;
}

i18n_locale_t *i18n_locale_init(const char *locale_name) {
  locale_t locale;
  i18n_locale_t *locale_struct;

  locale_struct = (i18n_locale_t *)sqlite3_malloc(sizeof(i18n_locale_t));
  if (locale_struct == NULL) {
    return NULL;
  }

  locale = newlocale(LC_ALL, "C", NULL);
  if (locale == NULL) {
    sqlite3_free(locale_struct);
    return NULL;
  }

  locale_struct->locale = locale;
  return locale_struct;
}

void i18n_locale_destroy(i18n_locale_t *locale) {
  if (locale != NULL) {
    freelocale(locale->locale);
    locale->locale = NULL;

    sqlite3_free(locale);
  }
}

#elif defined(LOCALE_USE_SET_LOCALE)

#include <locale.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct i18n_locale {
  int dummy;
};

static struct i18n_locale DUMMY_LOCALE;

/*
 * Fallback for platforms that offer neither _create_locale nor newlocale. setlocale changes
 * process-global state, so this cannot be made thread safe; parsing WKT from several threads at
 * once races here. The platforms that reach this branch have no per-thread alternative.
 */
double i18n_strtod(const char *nptr, char **endptr, i18n_locale_t *locale) {
  /* setlocale returns the locale that is now in effect, so the current one has to be queried
     before switching. Passing its return value to the restore call put "C" back instead of what
     the host had chosen, permanently changing the process's numeric locale. The name also has to
     be copied: setlocale returns static storage that the next call may overwrite. */
  char saved[64];
  const char *current = setlocale(LC_NUMERIC, NULL);
  int restore = 0;

  if (current != NULL && strlen(current) < sizeof(saved)) {
    strcpy(saved, current);
    restore = 1;
  }

  setlocale(LC_NUMERIC, "C");
  double result = strtod(nptr, endptr);

  if (restore) {
    setlocale(LC_NUMERIC, saved);
  }

  return result;
}

/* Same process-global setlocale dance, and therefore the same documented race, as i18n_strtod
   above. */
int i18n_snprintf(char *out, size_t length, i18n_locale_t *locale, const char *format, ...) {
  char saved[64];
  const char *current = setlocale(LC_NUMERIC, NULL);
  int restore = 0;

  if (current != NULL && strlen(current) < sizeof(saved)) {
    strcpy(saved, current);
    restore = 1;
  }

  setlocale(LC_NUMERIC, "C");
  va_list args;
  va_start(args, format);
  int result = vsnprintf(out, length, format, args);
  va_end(args);

  if (restore) {
    setlocale(LC_NUMERIC, saved);
  }

  return result;
}

i18n_locale_t *i18n_locale_init(const char *locale_name) { return &DUMMY_LOCALE; }

void i18n_locale_destroy(i18n_locale_t *locale) {}

#else

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>

struct i18n_locale {
  int dummy;
};

static struct i18n_locale DUMMY_LOCALE;

double i18n_strtod(const char *nptr, char **endptr, i18n_locale_t *locale) { return strtod(nptr, endptr); }

int i18n_snprintf(char *out, size_t length, i18n_locale_t *locale, const char *format, ...) {
  va_list args;
  va_start(args, format);
  int result = vsnprintf(out, length, format, args);
  va_end(args);
  return result;
}

i18n_locale_t *i18n_locale_init(const char *locale_name) { return &DUMMY_LOCALE; }

void i18n_locale_destroy(i18n_locale_t *locale) {}

#endif
