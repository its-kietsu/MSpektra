/*
 * MSpektra.exe - launcher of the portable MSpektra package.
 * Starts  python\pythonw.exe -s _portable\launch_unidec.py [arguments]
 * from the folder this exe is in, without a console window.
 * Deconvolution engine UniDec: Marty et al., Anal. Chem. 2015, DOI 10.1021/acs.analchem.5b00140
 */
#define WIN32_LEAN_AND_MEAN


#include <windows.h>
#include <shellapi.h>
#include <wchar.h>

#define CMD_MAX 32768

static void show_error(const wchar_t *msg, const wchar_t *detail)
{
    wchar_t text[4096];
    _snwprintf(text, 4096, L"%ls\n\n%ls", msg, detail ? detail : L"");
    text[4095] = 0;
    MessageBoxW(NULL, text, L"MSpektra", MB_OK | MB_ICONERROR);
}

/* Skip argv[0] in the raw command line (handles quoted paths). */
static const wchar_t *skip_program_name(const wchar_t *cmd)
{
    if (*cmd == L'"') {
        cmd++;
        while (*cmd && *cmd != L'"') cmd++;
        if (*cmd == L'"') cmd++;
    } else {
        while (*cmd && *cmd != L' ' && *cmd != L'\t') cmd++;
    }
    while (*cmd == L' ' || *cmd == L'\t') cmd++;
    return cmd;
}

int WINAPI wWinMain(HINSTANCE hInst, HINSTANCE hPrev, PWSTR pCmdLine, int nShow)
{
    wchar_t dir[MAX_PATH * 4], pyw[MAX_PATH * 4], script[MAX_PATH * 4];
    static wchar_t cmd[CMD_MAX];
    DWORD n = GetModuleFileNameW(NULL, dir, MAX_PATH * 4);
    if (n == 0 || n >= MAX_PATH * 4) {
        show_error(L"Could not determine the program folder.", NULL);
        return 1;
    }
    wchar_t *slash = wcsrchr(dir, L'\\');
    if (slash) *slash = 0;

    _snwprintf(pyw, MAX_PATH * 4, L"%ls\\python\\pythonw.exe", dir);
    _snwprintf(script, MAX_PATH * 4, L"%ls\\_portable\\launch_unidec.py", dir);
    if (GetFileAttributesW(pyw) == INVALID_FILE_ATTRIBUTES ||
        GetFileAttributesW(script) == INVALID_FILE_ATTRIBUTES) {
        show_error(L"MSpektra.exe must stay in its program folder,\n"
                   L"next to the 'python' and '_portable' folders.\n\nMissing:", pyw);
        return 1;
    }

    const wchar_t *args = skip_program_name(GetCommandLineW());
    int len = _snwprintf(cmd, CMD_MAX, L"\"%ls\" -s \"%ls\"%ls%ls",
                         pyw, script, (*args ? L" " : L""), args);
    if (len < 0 || len >= CMD_MAX) {
        show_error(L"Command line too long.", NULL);
        return 1;
    }
    cmd[CMD_MAX - 1] = 0;

    STARTUPINFOW si;
    PROCESS_INFORMATION pi;
    ZeroMemory(&si, sizeof(si));
    si.cb = sizeof(si);
    ZeroMemory(&pi, sizeof(pi));
    if (!CreateProcessW(pyw, cmd, NULL, NULL, FALSE, 0, NULL, dir, &si, &pi)) {
        DWORD err = GetLastError();
        wchar_t buf[512];
        FormatMessageW(FORMAT_MESSAGE_FROM_SYSTEM | FORMAT_MESSAGE_IGNORE_INSERTS, NULL, err, 0,
                       buf, 512, NULL);
        show_error(L"MSpektra could not be started (it may be blocked by security software).\n"
                   L"Try MSpektra.bat or MSpektra (console).bat instead.", buf);
        return 1;
    }
    AllowSetForegroundWindow(pi.dwProcessId);
    CloseHandle(pi.hThread);
    CloseHandle(pi.hProcess);
    return 0;
}
