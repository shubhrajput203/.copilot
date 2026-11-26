/*
 * LifeLine - Blood donor registration and urgent requests (console app in C)
 *
 * Build:  gcc -std=c99 -Wall -O2 blood_donor.c -o blood_donor -lm
 * Run:    ./blood_donor
 *
 * Data is saved in donors.dat and requests.dat in the current folder.
 */
#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <time.h>

#define MAX_DONORS 1000
#define MAX_REQS   1000
#define NAME_LEN   48
#define AREA_LEN   64
#define PHONE_LEN  20
#define NOTE_LEN   96
#define NO_COORD   999.0
#define WAIT_DAYS  90

static const char *GROUPS[8] = {"A+", "A-", "B+", "B-", "AB+", "AB-", "O+", "O-"};
static const char *AVAIL[3]  = {"Available now", "Weekends only", "Unavailable"};
static const char *URGENCY[3] = {"Critical", "Urgent", "Normal"};
static const char *STATUS[4] = {"Open", "In progress", "Fulfilled", "Cancelled"};

typedef struct {
    int id;
    char name[NAME_LEN], phone[PHONE_LEN], area[AREA_LEN];
    int group, avail;          /* group: index into GROUPS, avail: index into AVAIL */
    int days_since_donation;   /* -1 = never or unknown */
    double lat, lng;           /* NO_COORD when not given */
} Donor;

typedef struct {
    int id;                    /* ticket number, shown as BR-<id> */
    char patient[NAME_LEN], phone[PHONE_LEN], area[AREA_LEN], notes[NOTE_LEN];
    int group, units, urgency, status;
    double lat, lng;
    time_t created, updated;
} Request;

static Donor donors[MAX_DONORS];
static Request reqs[MAX_REQS];
static int n_donors = 0, n_reqs = 0;

/* ---------- input helpers ---------- */

static void read_line(const char *prompt, char *buf, size_t size) {
    printf("%s", prompt);
    fflush(stdout);
    if (!fgets(buf, (int)size, stdin)) { buf[0] = '\0'; return; }
    buf[strcspn(buf, "\r\n")] = '\0';
    for (char *p = buf; *p; p++) if (*p == '|') *p = '/';
}

static int read_int(const char *prompt, int lo, int hi) {
    char buf[64];
    for (;;) {
        read_line(prompt, buf, sizeof buf);
        char *end;
        long v = strtol(buf, &end, 10);
        if (end != buf && *end == '\0' && v >= lo && v <= hi) return (int)v;
        printf("  Enter a number from %d to %d.\n", lo, hi);
    }
}

static void read_required(const char *prompt, char *buf, size_t size) {
    do {
        read_line(prompt, buf, size);
        if (!buf[0]) printf("  This field is required.\n");
    } while (!buf[0]);
}

static int pick(const char *title, const char **opts, int n) {
    printf("%s\n", title);
    for (int i = 0; i < n; i++) printf("  %d) %s\n", i + 1, opts[i]);
    return read_int("Choose: ", 1, n) - 1;
}

/* Reads a latitude/longitude pair. Blank = skip, so area text is used instead. */
static void read_location(double *lat, double *lng) {
    char buf[96];
    *lat = *lng = NO_COORD;
    read_line("Location as \"latitude longitude\" (e.g. 12.9716 77.5946), or press Enter to skip: ", buf, sizeof buf);
    if (!buf[0]) return;
    double a, b;
    if (sscanf(buf, "%lf %lf", &a, &b) == 2 && fabs(a) <= 90 && fabs(b) <= 180) {
        *lat = a; *lng = b;
    } else {
        printf("  Could not read that location. Using area text only.\n");
    }
}

/* ---------- compatibility and distance ---------- */

/* Group index -> ABO bits (A=1, B=2, AB=3, O=0) and Rh (1 = positive). */
static int abo_of(int g) { static const int t[8] = {1, 1, 2, 2, 3, 3, 0, 0}; return t[g]; }
static int rh_of(int g)  { return g % 2 == 0; }

/* A donor can give to a patient when the donor's antigens are a subset of the
 * patient's and a negative patient only receives Rh negative blood. */
static int can_donate(int donor_group, int patient_group) {
    int d = abo_of(donor_group), p = abo_of(patient_group);
    if ((d & ~p) != 0) return 0;
    if (rh_of(donor_group) && !rh_of(patient_group)) return 0;
    return 1;
}

static double haversine_km(double la1, double lo1, double la2, double lo2) {
    const double R = 6371.0, r = 3.14159265358979323846 / 180.0;
    double dla = (la2 - la1) * r, dlo = (lo2 - lo1) * r;
    double h = sin(dla / 2) * sin(dla / 2) +
               cos(la1 * r) * cos(la2 * r) * sin(dlo / 2) * sin(dlo / 2);
    return 2 * R * asin(sqrt(h));
}

static int has_coord(double lat) { return lat != NO_COORD; }

static int ci_contains(const char *hay, const char *needle) {
    if (!*needle) return 0;
    size_t n = strlen(needle);
    for (; *hay; hay++) {
        size_t i = 0;
        while (i < n && hay[i] && tolower((unsigned char)hay[i]) == tolower((unsigned char)needle[i])) i++;
        if (i == n) return 1;
    }
    return 0;
}

static int donor_eligible(const Donor *d) {
    return d->days_since_donation < 0 || d->days_since_donation >= WAIT_DAYS;
}

/* ---------- matching (compatible groups, nearest first) ---------- */

typedef struct { const Donor *d; double km; int same_area; } Match;

static int cmp_match(const void *a, const void *b) {
    const Match *x = a, *y = b;
    int xc = has_coord(x->km), yc = has_coord(y->km);
    if (xc && yc) return (x->km > y->km) - (x->km < y->km);
    if (xc != yc) return xc ? -1 : 1;
    if (x->same_area != y->same_area) return y->same_area - x->same_area;
    return 0;
}

static int find_matches(const Request *r, Match *out) {
    int n = 0;
    for (int i = 0; i < n_donors; i++) {
        const Donor *d = &donors[i];
        if (!can_donate(d->group, r->group)) continue;
        if (d->avail == 2 || !donor_eligible(d)) continue;
        out[n].d = d;
        out[n].km = (has_coord(d->lat) && has_coord(r->lat))
                        ? haversine_km(d->lat, d->lng, r->lat, r->lng) : NO_COORD;
        out[n].same_area = ci_contains(r->area, d->area) || ci_contains(d->area, r->area);
        n++;
    }
    qsort(out, n, sizeof *out, cmp_match);
    return n;
}

/* ---------- storage ---------- */

static void save_all(void) {
    FILE *f = fopen("donors.dat", "wb");
    if (f) { fwrite(&n_donors, sizeof n_donors, 1, f); fwrite(donors, sizeof(Donor), n_donors, f); fclose(f); }
    else perror("Could not save donors.dat");
    f = fopen("requests.dat", "wb");
    if (f) { fwrite(&n_reqs, sizeof n_reqs, 1, f); fwrite(reqs, sizeof(Request), n_reqs, f); fclose(f); }
    else perror("Could not save requests.dat");
}

static void load_all(void) {
    FILE *f = fopen("donors.dat", "rb");
    if (f) {
        if (fread(&n_donors, sizeof n_donors, 1, f) != 1 || n_donors < 0 || n_donors > MAX_DONORS) n_donors = 0;
        else if ((int)fread(donors, sizeof(Donor), n_donors, f) != n_donors) n_donors = 0;
        fclose(f);
    }
    f = fopen("requests.dat", "rb");
    if (f) {
        if (fread(&n_reqs, sizeof n_reqs, 1, f) != 1 || n_reqs < 0 || n_reqs > MAX_REQS) n_reqs = 0;
        else if ((int)fread(reqs, sizeof(Request), n_reqs, f) != n_reqs) n_reqs = 0;
        fclose(f);
    }
}

static int next_donor_id(void) { int m = 0; for (int i = 0; i < n_donors; i++) if (donors[i].id > m) m = donors[i].id; return m + 1; }
static int next_ticket_id(void) { int m = 1000; for (int i = 0; i < n_reqs; i++) if (reqs[i].id > m) m = reqs[i].id; return m + 1; }

/* ---------- actions ---------- */

static void register_donor(void) {
    if (n_donors >= MAX_DONORS) { printf("Donor list is full.\n"); return; }
    Donor *d = &donors[n_donors];
    memset(d, 0, sizeof *d);
    d->id = next_donor_id();
    printf("\n== Register as a donor ==\n");
    read_required("Full name: ", d->name, sizeof d->name);
    read_required("Phone number: ", d->phone, sizeof d->phone);
    d->group = pick("Blood group:", GROUPS, 8);
    d->avail = pick("Availability:", AVAIL, 3);
    read_required("Your area (neighbourhood or locality): ", d->area, sizeof d->area);
    d->days_since_donation = read_int("Days since last donation (0 to 3650, or -1 if never): ", -1, 3650);
    read_location(&d->lat, &d->lng);
    n_donors++;
    save_all();
    printf("Registered %s (%s) as donor #%d.\n", d->name, GROUPS[d->group], d->id);
}

static void raise_request(void) {
    if (n_reqs >= MAX_REQS) { printf("Request list is full.\n"); return; }
    Request *r = &reqs[n_reqs];
    memset(r, 0, sizeof *r);
    r->id = next_ticket_id();
    printf("\n== Raise an urgent request ==\n");
    read_required("Patient name: ", r->patient, sizeof r->patient);
    r->group = pick("Blood group needed:", GROUPS, 8);
    r->units = read_int("Units needed (1 to 20): ", 1, 20);
    r->urgency = pick("Urgency:", URGENCY, 3);
    read_required("Hospital and area: ", r->area, sizeof r->area);
    read_required("Contact number: ", r->phone, sizeof r->phone);
    read_line("Notes (ward, doctor, time needed by): ", r->notes, sizeof r->notes);
    read_location(&r->lat, &r->lng);
    r->status = 0;
    r->created = r->updated = time(NULL);
    n_reqs++;
    save_all();
    printf("Ticket BR-%d raised. Use option 4 to see matching donors.\n", r->id);
}

static void fmt_time(time_t t, char *buf, size_t n) {
    struct tm *tm = localtime(&t);
    strftime(buf, n, "%d %b %H:%M", tm);
}

static void print_request(const Request *r) {
    char c[32], u[32];
    fmt_time(r->created, c, sizeof c);
    fmt_time(r->updated, u, sizeof u);
    printf("BR-%d | %-11s | %-8s | needs %-3s x%d | %s\n", r->id, STATUS[r->status], URGENCY[r->urgency],
           GROUPS[r->group], r->units, r->patient);
    printf("        %s | %s | raised %s, updated %s\n", r->area, r->phone, c, u);
    if (r->notes[0]) printf("        Notes: %s\n", r->notes);
}

static int by_urgency_then_new(const void *a, const void *b) {
    const Request *x = *(Request *const *)a, *y = *(Request *const *)b;
    if (x->urgency != y->urgency) return x->urgency - y->urgency;
    return y->id - x->id;
}

static void list_requests(void) {
    printf("\nShow which tickets?\n  1) Active (Open and In progress)\n  2) All\n  3) One status\n");
    int mode = read_int("Choose: ", 1, 3), only = -1;
    if (mode == 3) only = pick("Status:", STATUS, 4);
    Request *view[MAX_REQS];
    int n = 0;
    for (int i = 0; i < n_reqs; i++) {
        int st = reqs[i].status;
        if ((mode == 1 && st > 1) || (mode == 3 && st != only)) continue;
        view[n++] = &reqs[i];
    }
    if (!n) { printf("No requests here yet. Choose \"Raise urgent request\" to open the first ticket.\n"); return; }
    qsort(view, n, sizeof *view, by_urgency_then_new);
    printf("\n");
    for (int i = 0; i < n; i++) print_request(view[i]);
}

static Request *find_request(int id) {
    for (int i = 0; i < n_reqs; i++) if (reqs[i].id == id) return &reqs[i];
    return NULL;
}

static Request *ask_request(void) {
    if (!n_reqs) { printf("No requests yet.\n"); return NULL; }
    int id = read_int("Ticket number (the digits after BR-): ", 1, 999999);
    Request *r = find_request(id);
    if (!r) printf("No ticket BR-%d found.\n", id);
    return r;
}

static void show_matches(void) {
    Request *r = ask_request();
    if (!r) return;
    static Match m[MAX_DONORS];
    int n = find_matches(r, m);
    printf("\n");
    print_request(r);
    printf("\nAccepted donor groups:");
    for (int g = 0; g < 8; g++) if (can_donate(g, r->group)) printf(" %s", GROUPS[g]);
    printf("\n%d compatible, available donor%s, nearest first:\n", n, n == 1 ? "" : "s");
    if (!n) { printf("  None yet. Share the ticket or ask donors to register.\n"); return; }
    for (int i = 0; i < n; i++) {
        char dist[32];
        if (has_coord(m[i].km)) {
            if (m[i].km < 1) snprintf(dist, sizeof dist, "%d m", (int)(m[i].km * 1000));
            else snprintf(dist, sizeof dist, "%.1f km", m[i].km);
        } else snprintf(dist, sizeof dist, "%s", m[i].same_area ? "same area" : "unknown");
        printf("  %2d. %-8s %-3s %-22s %-14s %s | %s\n", i + 1, dist, GROUPS[m[i].d->group], m[i].d->name,
               m[i].d->phone, m[i].d->area, AVAIL[m[i].d->avail]);
    }
}

static void update_status(void) {
    Request *r = ask_request();
    if (!r) return;
    printf("\nCurrent status of BR-%d: %s\n", r->id, STATUS[r->status]);
    int s = pick("New status:", STATUS, 4);
    if (s == r->status) { printf("Status is already %s.\n", STATUS[s]); return; }
    r->status = s;
    r->updated = time(NULL);
    save_all();
    printf("Ticket BR-%d is now %s.\n", r->id, STATUS[s]);
}

static void list_donors(void) {
    if (!n_donors) { printf("No donors yet. Register the first donor.\n"); return; }
    printf("\n%d registered donor%s (eligible again %d days after a donation)\n", n_donors, n_donors == 1 ? "" : "s", WAIT_DAYS);
    for (int i = 0; i < n_donors; i++) {
        const Donor *d = &donors[i];
        printf("#%-3d %-3s %-22s %-14s %-18s %s", d->id, GROUPS[d->group], d->name, d->phone, d->area, AVAIL[d->avail]);
        if (!donor_eligible(d)) printf(" (not eligible for %d more days)", WAIT_DAYS - d->days_since_donation);
        printf("\n");
    }
}

static void update_availability(void) {
    if (!n_donors) { printf("No donors yet.\n"); return; }
    int id = read_int("Donor number: ", 1, 999999);
    for (int i = 0; i < n_donors; i++) {
        if (donors[i].id != id) continue;
        printf("%s is currently: %s\n", donors[i].name, AVAIL[donors[i].avail]);
        donors[i].avail = pick("New availability:", AVAIL, 3);
        save_all();
        printf("Availability updated.\n");
        return;
    }
    printf("No donor #%d found.\n", id);
}

/* ---------- main ---------- */

int main(void) {
    load_all();
    for (;;) {
        printf("\n=========== LifeLine ===========\n"
               "1) Register as donor\n"
               "2) Raise urgent request\n"
               "3) View requests\n"
               "4) Find compatible donors for a request\n"
               "5) Update ticket status\n"
               "6) List donors\n"
               "7) Update donor availability\n"
               "0) Exit\n");
        switch (read_int("Choose: ", 0, 7)) {
            case 1: register_donor(); break;
            case 2: raise_request(); break;
            case 3: list_requests(); break;
            case 4: show_matches(); break;
            case 5: update_status(); break;
            case 6: list_donors(); break;
            case 7: update_availability(); break;
            case 0: save_all(); printf("Saved. Stay safe.\n"); return 0;
        }
    }
}