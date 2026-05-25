self.Module = self.Module || {};
self.Module.preRun = self.Module.preRun || [];

self.Module.preRun.push(function() {
    addRunDependency('cocoa_fs_init');

    async function initFS() {
        try {
            const mapRes = await fetch('cocoa-fs.json');
            const fileList = await mapRes.json();
            const vfsBase = '/src/CoCoA-5';

            var createdDirs = {};

            fileList.forEach(function(appPath) {
                var parts = appPath.split('/');
                var fileName = parts.pop();
                var dirPath = parts.join('/');
                
                var vfsDir = vfsBase + (dirPath ? '/' + dirPath : '');

                if (!createdDirs[vfsDir]) {
                    FS.mkdirTree(vfsDir);
                    createdDirs[vfsDir] = true;
                }

                var fetchPath = appPath; 
                FS.createLazyFile(vfsDir, fileName, fetchPath, true, false);
            });

            removeRunDependency('cocoa_fs_init');
        } catch(err) {
            console.error("Failed to initialize CoCoA FS:", err);
            removeRunDependency('cocoa_fs_init');
        }
    }

    initFS();
});