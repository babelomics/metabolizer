#!/usr/bin/env node

'use strict';

const Vulcanize = require('vulcanize');
const async = require('async');
const path = require('path');
const shell = require('shelljs');
const fs = require('fs');

var name = "metabolizer";
var element = "metabolizer-element";

// Default
var bp = path.join(__dirname, "build");
var indexHTML = path.join(__dirname, name + '-index.html');
var buildIndexHTML = path.join(bp, 'index.html');
var elementHTML = path.join(__dirname, element + '.html');
var buildElementHTML = path.join(bp, element + '.html');

var vulcan = new Vulcanize({
    inlineScripts: true,
    inlineCss: true,
    stripComments: true
});

async.waterfall([
    function (cb) {
        shell.rm('-rf', bp);
        shell.mkdir('-p', bp);
        shell.mkdir('-p', path.join(bp, "fontawesome"));
        shell.mkdir('-p', path.join(bp, "webcomponentsjs"));

        cb(null);
    },
    function (cb) {
        vulcan.process(elementHTML, function (err, inlinedHtml) {
            fs.writeFile(buildElementHTML, inlinedHtml, function (err) {
                if (err) {
                    cb(err);
                } else {
                    cb(null);
                }
            });
        });
    },
    function (cb) {
        shell.cp('-r', indexHTML, buildIndexHTML);
        shell.cp('-r', path.join(__dirname, 'conf/'), bp);
        shell.cp('-r', path.join(__dirname, 'images/'), bp);
        shell.cp('-r', path.join(__dirname, 'bower_components', 'stevia-elements', 'fonts'), bp);
        shell.cp('-r', path.join(__dirname, 'bower_components', 'stevia-elements', 'css'), bp);
        shell.cp('-r', path.join(__dirname, 'bower_components', 'fontawesome', 'css'), path.join(bp, "fontawesome/"));
        shell.cp('-r', path.join(__dirname, 'bower_components', 'fontawesome', 'fonts'), path.join(bp, "fontawesome/"));
        shell.cp('-r', path.join(__dirname, 'bower_components', 'webcomponentsjs', '*.min.js'), path.join(bp, "webcomponentsjs/"));

        // fix index.html paths
        shell.sed('-i', 'bower_components/stevia-elements/', '', buildIndexHTML);
        shell.sed('-i', 'bower_components/', '', buildIndexHTML);

        cb(null);
    }
], function (err) {
    if (err) {
        console.log(err);
    }
    console.log("Done.");
});
