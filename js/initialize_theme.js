(function () {
    const currentTheme = localStorage.getItem('theme');
    if (currentTheme) {
        document.documentElement.setAttribute('data-theme', currentTheme);
    } else {
        // NOTE: disabled, should be always light unless switched,
        // because of images that are not suitable for dark mode
        // const isSystemDark = window.matchMedia('(prefers-color-scheme: dark)').matches;
        // document.documentElement.setAttribute('data-theme', isSystemDark ? 'dark' : 'light');

        document.documentElement.setAttribute('data-theme', 'light');
    }
})();
